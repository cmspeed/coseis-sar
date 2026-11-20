"""
Lossless repack of a COSEIS_SAR netCDF for the archive, checked against its original.

The layers that DockerizedTopsApp's add_2d_layer writes uncompressed and contiguous (ionosphere,
ionosphere burst ramps, unfiltered coherence, range and azimuth offsets, offset SNR) get the
chunking of the standard layers (losDisplacement's), shuffle and gzip 4, as the standard layers
have. Everything else is kept exactly as it is (rebuild()): the new file is created with the
original's file-creation properties, every other object is copied with H5Ocopy (object header,
attributes in their stored order, data chunks as they are), and only the compressed layers (and
the groups above them) are created anew. h5repack cannot do this: it drops the groups' attribute
creation-order flag, which changes the order netCDF lists attributes in. A product with nothing
left to compress (after the add_2d_layer fix) is copied the same way.

The file is written to a hidden .part name, flushed and dropped from the page cache (so the check
reads it back from disk), checked by check_pair(), and only then renamed. Other files (the PNG)
are copied with direct I/O and checked by SHA-256. The source is only read.

This is the code that built COSEIS-ONE_STOP_SHOP_COMPRESSED (2026-10-28), run here one file at a
time. It needs h5py, numpy and ncdump (the topsApp environment), so the aria_coseis package never
imports it: aria_coseis.sar.products runs it as a script.

    python -I repack.py SRC DST

Prints one JSON record; exit status 0 if the file was placed and passed every check.
"""

from __future__ import annotations

import difflib
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

import h5py
import numpy as np

BIN = Path(sys.executable).parent
ADDITIONAL_LAYERS = (
    "science/grids/corrections/derived/ionosphere/ionosphere",
    "science/grids/corrections/derived/ionosphereBurstRamps/ionosphereBurstRamps",
    "science/grids/data/unfilteredCoherence",
    "science/grids/data/azimuthPixelOffsets",
    "science/grids/data/rangePixelOffsets",
    "science/grids/data/denseOffsetsSNR",
)
REFERENCE_LAYER = "science/grids/data/losDisplacement"
NC_STORAGE_ATTRS = re.compile(
    r'^(_Storage = "chunked"|_ChunkSizes = \d+(, \d+)*|_Shuffle = "true"|_DeflateLevel = 4) ;$'
)


def targets(path: Path) -> dict[str, tuple[int, ...]]:
    """Layers to compress, with their new chunk shape (the reference layer's, clipped to the shape)."""
    with h5py.File(path, "r") as f:
        ref = f[REFERENCE_LAYER].chunks
        found: dict[str, tuple[int, ...]] = {}
        for name in ADDITIONAL_LAYERS:
            obj = f.get(name)
            if (
                isinstance(obj, h5py.Dataset)
                and obj.ndim == 2
                and obj.chunks is None
                and obj.compression is None
            ):
                found[name] = tuple(min(c, s) for c, s in zip(ref, obj.shape))
    return found


def plain(f: h5py.File, value: Any) -> Any:
    """An attribute value with object references replaced by the names they point to."""
    if isinstance(value, h5py.Reference):
        # Null, dangling and unlinked-object references compare equal: the originals' matchup scale
        # has two in its REFERENCE_LIST to variables that packaging deleted (unlinked objects, which
        # h5repack does not copy, so in the copy they point nowhere)
        if not value:
            return ("ref", "<unusable>")
        try:
            return ("ref", f[value].name or "<unusable>")
        except (KeyError, ValueError):
            return ("ref", "<unusable>")
    if isinstance(value, np.ndarray):
        if value.dtype.names:  # e.g. REFERENCE_LIST: (dataset reference, dimension index)
            return [tuple(plain(f, v) for v in row) for row in value.tolist()]
        if value.dtype.kind == "O":
            return [plain(f, v) for v in value.ravel().tolist()] + [value.shape]
        return ("array", value.dtype.str, value.shape, value.tobytes())
    if isinstance(value, (list, tuple)):
        return [plain(f, v) for v in value]
    if isinstance(value, np.generic):
        return ("scalar", value.dtype.str, value.tobytes())
    return value


HIDDEN = (
    "DIMENSION_LIST",
    "REFERENCE_LIST",
)  # HDF5 dimension-scale bookkeeping (not netCDF attributes)


def stored_order(oid: Any) -> list[str]:
    """Attribute names in the order netCDF lists them: creation order where it is tracked, otherwise
    the order they are stored in (message order when compact, name-hash order when dense)."""
    tracked = oid.get_create_plist().get_attr_creation_order() & h5py.h5p.CRT_ORDER_TRACKED
    names: list[str] = []
    h5py.h5a.iterate(
        oid,
        lambda name, *_: names.append(name.decode()),
        index_type=h5py.h5.INDEX_CRT_ORDER if tracked else h5py.h5.INDEX_NAME,
        order=h5py.h5.ITER_INC if tracked else h5py.h5.ITER_NATIVE,
    )
    return [n for n in names if n not in HIDDEN]


def object_ids(f: h5py.File) -> dict[str, Any]:
    ids = {"": h5py.h5g.open(f.id, b"/")}

    def add(name: str, obj: Any) -> None:  # must return None, or visititems stops
        ids[name] = obj.id

    f.visititems(add)
    return ids


def copy_attr(src: Any, dst_loc: Any, name: bytes, refs_as_null: bool = False) -> None:
    """Create attribute name on dst_loc with src's (an h5py AttrID) type, shape and bytes: a raw copy
    in its own file type, except variable-length data and references, which go through h5py. With
    refs_as_null, references are written null (rebuild() sets them afterwards, in place)."""
    tid, space = src.get_type(), src.get_space()
    new = h5py.h5a.create(dst_loc, name, tid, space)
    if (
        tid.detect_class(h5py.h5t.REFERENCE)
        or tid.detect_class(h5py.h5t.VLEN)
        or (tid.get_class() == h5py.h5t.STRING and tid.is_variable_str())
    ):
        buf = np.empty(src.shape, dtype=src.dtype)
        src.read(buf)
        if refs_as_null:
            buf = map_refs(buf, lambda r: h5py.Reference())
        new.write(buf)
    else:
        buf = np.empty(src.shape, dtype=np.dtype(f"V{tid.get_size()}"))
        src.read(buf, mtype=tid)
        new.write(buf, mtype=tid)


def map_refs(value: Any, fn) -> Any:
    """value (as h5py reads an attribute) with every object reference r replaced by fn(r)."""
    if isinstance(value, h5py.Reference):
        return fn(value)
    if isinstance(value, np.ndarray):
        out = np.empty(value.shape, dtype=value.dtype)
        if value.dtype.names:
            for field in value.dtype.names:
                out[field] = map_refs(value[field], fn)
        elif value.dtype.kind == "O":
            for i, v in np.ndenumerate(value):
                out[i] = map_refs(v, fn)
        else:
            out[...] = value
        return out
    return value


def attr_names(oid: Any) -> list[str]:
    """All attribute names (hidden ones too) in stored order (see stored_order)."""
    tracked = oid.get_create_plist().get_attr_creation_order() & h5py.h5p.CRT_ORDER_TRACKED
    names: list[str] = []
    h5py.h5a.iterate(
        oid,
        lambda name, *_: names.append(name.decode()),
        index_type=h5py.h5.INDEX_CRT_ORDER if tracked else h5py.h5.INDEX_NAME,
        order=h5py.h5.ITER_INC if tracked else h5py.h5.ITER_NATIVE,
    )
    return names


def links(gid: Any) -> list[bytes]:
    """A group's link names in creation order where it is tracked (else by name)."""
    tracked = gid.get_create_plist().get_link_creation_order() & h5py.h5p.CRT_ORDER_TRACKED
    names: list[bytes] = []
    gid.links.iterate(
        lambda name: names.append(name),
        idx_type=h5py.h5.INDEX_CRT_ORDER if tracked else h5py.h5.INDEX_NAME,
        order=h5py.h5.ITER_INC,
    )
    return names


# Property lists are made fresh with the original's settings: the ones HDF5 returns for an existing
# object carry that file's internals (a group's link storage), and reusing them corrupts the copy.
def fresh_fcpl(f: h5py.File) -> Any:
    src = f.id.get_create_plist()
    new = h5py.h5p.create(h5py.h5p.FILE_CREATE)
    default = h5py.h5p.create(h5py.h5p.FILE_CREATE)
    for get in (
        "get_sizes",
        "get_userblock",
        "get_file_space_strategy",
        "get_file_space_page_size",
    ):
        if getattr(src, get)() != getattr(default, get)():
            raise RuntimeError(
                f"file creation property {get[4:]} is not the default: {getattr(src, get)()}"
            )
    root = h5py.h5g.open(f.id, b"/").get_create_plist()  # the root group's flags are reported here
    for prop in ("link_creation_order", "attr_creation_order", "obj_track_times"):
        getattr(new, "set_" + prop)(getattr(root, "get_" + prop)())
    new.set_attr_phase_change(*root.get_attr_phase_change())
    return new


def fresh_gcpl(gid: Any) -> Any:
    src = gid.get_create_plist()
    new = h5py.h5p.create(h5py.h5p.GROUP_CREATE)
    for prop in ("link_creation_order", "attr_creation_order", "obj_track_times"):
        getattr(new, "set_" + prop)(getattr(src, "get_" + prop)())
    new.set_attr_phase_change(*src.get_attr_phase_change())
    return new


def fresh_dcpl(dsid: Any, chunks: tuple[int, ...] | None) -> Any:
    """The original's dataset settings; with chunks, chunked with shuffle and deflate 4 instead of
    the original's layout and filters."""
    src = dsid.get_create_plist()
    if src.get_external_count():
        raise RuntimeError("dataset with external storage")
    new = h5py.h5p.create(h5py.h5p.DATASET_CREATE)
    if chunks:
        if src.get_nfilters():
            raise RuntimeError("layer to compress already has filters")
        new.set_chunk(chunks)
        new.set_shuffle()
        new.set_deflate(4)
        new.set_alloc_time(h5py.h5d.ALLOC_TIME_INCR)
    else:
        layout = src.get_layout()
        if layout == h5py.h5d.CHUNKED:
            new.set_chunk(src.get_chunk())
        else:
            new.set_layout(layout)
        for i in range(src.get_nfilters()):
            code, flags, values, _ = src.get_filter(i)
            new.set_filter(code, flags, values)
        new.set_alloc_time(src.get_alloc_time())
    new.set_fill_time(src.get_fill_time())
    if src.fill_value_defined() == h5py.h5d.FILL_VALUE_USER_DEFINED:
        fill = np.zeros(
            (1,), dtype=h5py.Dataset(dsid).dtype
        )  # h5py needs one element for vlen types
        src.get_fill_value(fill)
        new.set_fill_value(fill)
    for prop in ("attr_creation_order", "obj_track_times"):
        getattr(new, "set_" + prop)(getattr(src, "get_" + prop)())
    new.set_attr_phase_change(*src.get_attr_phase_change())
    return new


def rebuild(src: Path, dst: Path, tg: dict[str, tuple[int, ...]]) -> None:
    """Write dst: src with the layers tg compressed and everything else as it is (module docstring)."""
    with h5py.File(src, "r") as fa:
        fcpl = fresh_fcpl(fa)
        fapl = h5py.h5p.create(h5py.h5p.FILE_ACCESS)
        fapl.set_libver_bounds(*fa.id.get_access_plist().get_libver_bounds())
        fid = h5py.h5f.create(str(dst).encode(), h5py.h5f.ACC_TRUNC, fcpl=fcpl, fapl=fapl)
        with h5py.File(fid) as fb:
            # Created anew: the layers to compress, and datasets whose fill value is variable-length
            # (the radarMetaData strings): H5Ocopy carries such a fill value over still pointing into
            # the source file's global heap, which breaks the copy
            special = dict(tg)

            def vlen_fill(name: str, obj: Any) -> None:
                if isinstance(obj, h5py.Dataset) and name not in special:
                    t = obj.id.get_type()
                    vlen = t.detect_class(h5py.h5t.VLEN) or (
                        t.get_class() == h5py.h5t.STRING and t.is_variable_str()
                    )
                    if (
                        vlen
                        and obj.id.get_create_plist().fill_value_defined()
                        == h5py.h5d.FILL_VALUE_USER_DEFINED
                    ):
                        special[name] = None

            fa.visititems(vlen_fill)
            above = {""} | {
                "/".join(t.split("/")[:i]) for t in special for i in range(1, t.count("/") + 1)
            }

            def copy_attrs(sid: Any, did: Any) -> None:
                for name in attr_names(sid):
                    copy_attr(
                        h5py.h5a.open(sid, name.encode()), did, name.encode(), refs_as_null=True
                    )

            def group(path: str, sgid: Any, dgid: Any) -> None:
                for name in links(sgid):
                    info = sgid.links.get_info(name)
                    if info.type != h5py.h5l.TYPE_HARD:
                        raise RuntimeError(f"/{path}/{name.decode()}: not a hard link")
                    lcpl = h5py.h5p.create(h5py.h5p.LINK_CREATE)
                    lcpl.set_char_encoding(info.cset)
                    child = f"{path}/{name.decode()}".lstrip("/")
                    if h5py.h5o.get_info(sgid, name).rc != 1:
                        raise RuntimeError(f"/{child}: object with several links")
                    if child in above:
                        sub = h5py.h5g.open(sgid, name)
                        new = h5py.h5g.create(dgid, name, lcpl=lcpl, gcpl=fresh_gcpl(sub))
                        copy_attrs(sub, new)
                        group(child, sub, new)
                    elif child in special:
                        sds = h5py.h5d.open(sgid, name)
                        dds = h5py.h5d.create(
                            dgid,
                            name,
                            sds.get_type(),
                            sds.get_space(),
                            dcpl=fresh_dcpl(sds, special[child]),
                            lcpl=lcpl,
                        )
                        copy_attrs(sds, dds)
                        a, b = h5py.Dataset(sds), h5py.Dataset(dds)
                        if a.shape == ():
                            b[()] = a[()]
                        else:
                            step = special[child][0] if special[child] else a.shape[0]
                            for r in range(0, a.shape[0], step):
                                b[r : r + step] = a[r : r + step]
                    else:
                        h5py.h5o.copy(sgid, name, dgid, name, lcpl=lcpl)

            root_a, root_b = h5py.h5g.open(fa.id, b"/"), h5py.h5g.open(fb.id, b"/")
            copy_attrs(root_a, root_b)
            group("", root_a, root_b)

            # References (dimension scales): H5Ocopy keeps the source's addresses, and copy_attr wrote
            # nulls, so every reference attribute is written again, in place (the same message, so
            # nothing moves), pointing to the object of the same name; references to nothing usable
            # in the source (unlinked objects) stay null
            def target(ref: h5py.Reference) -> h5py.Reference:
                if not ref:
                    return h5py.Reference()
                try:
                    name = fa[ref].name
                except (KeyError, ValueError):
                    return h5py.Reference()
                return fb[name].ref if name and name in fb else h5py.Reference()

            ids_a, ids_b = object_ids(fa), object_ids(fb)
            for n, ida in ids_a.items():
                idb = ids_b[n]
                if isinstance(ida, h5py.h5t.TypeID):
                    continue
                for name in attr_names(ida):
                    aid = h5py.h5a.open(ida, name.encode())
                    if not aid.get_type().detect_class(h5py.h5t.REFERENCE):
                        continue
                    buf = np.empty(aid.shape, dtype=aid.dtype)
                    aid.read(buf)
                    h5py.h5a.open(idb, name.encode()).write(map_refs(buf, target))


def attrs_compare(fa: h5py.File, a: Any, fb: h5py.File, b: Any) -> str:
    """ "same": the same attribute names, HDF5 types (H5Tequal), shapes and values, in the same order
    as netCDF lists them; "order": all that except the order; "differ": anything else.

    "order" happens where an original stores attributes without creation order (the extra layers and
    some metadata variables of isce2_topsapp 1.0.2 and 1.0.4 products): their order is then the
    physical order in the object header, which h5repack's copy does not keep (HDF5 places attributes
    by free space), and nothing can restore it short of rewriting the file through netCDF (nccopy),
    which changes much more (string types, fill-value properties, storage of small variables)."""
    ka, kb = list(a.attrs.keys()), list(b.attrs.keys())
    if sorted(ka) != sorted(kb):
        return "differ"
    for k in ka:
        ia, ib = a.attrs.get_id(k), b.attrs.get_id(k)
        if ia.get_type() != ib.get_type() or ia.shape != ib.shape:
            return "differ"
        if plain(fa, a.attrs[k]) != plain(fb, b.attrs[k]):
            return "differ"
    oa = a.id if not isinstance(a, h5py.File) else h5py.h5g.open(a.id, b"/")
    ob = b.id if not isinstance(b, h5py.File) else h5py.h5g.open(b.id, b"/")
    return "same" if stored_order(oa) == stored_order(ob) else "order"


def same_values(a: h5py.Dataset, b: h5py.Dataset) -> bool:
    """Bit-for-bit equality, read in blocks of whole chunk rows (each chunk is decompressed once)."""
    if a.shape == () or a.size == 0:
        x, y = a[()], b[()]
        return (
            np.array_equal(np.asarray(x), np.asarray(y))
            if a.dtype.kind == "O"
            else np.asarray(x).tobytes() == np.asarray(y).tobytes()
        )
    step = (b.chunks or a.chunks or (1024,))[0] if a.ndim > 1 else a.shape[0]
    for r in range(0, a.shape[0], step):
        x, y = a[r : r + step], b[r : r + step]
        if a.dtype.kind in "fiub" and not a.dtype.names:
            u = f"u{a.dtype.itemsize}" if a.dtype.itemsize in (1, 2, 4, 8) else None
            same = np.array_equal(x.view(u), y.view(u)) if u else x.tobytes() == y.tobytes()
        else:
            same = np.array_equal(x, y)
        if not same:
            return False
    return True


def storage(d: h5py.Dataset) -> tuple:
    return (d.chunks, d.compression, d.compression_opts, d.shuffle, d.fletcher32, d.scaleoffset)


def ncdump(path: Path, flags: str) -> list[str]:
    """ncdump output without its first line ("netcdf <file name> {"; the file being checked is .part)."""
    out = subprocess.run(
        [str(BIN / "ncdump"), flags, str(path)], capture_output=True, text=True, check=True
    )
    return out.stdout.splitlines()[1:]


def check_pair(orig: Path, new: Path, tg: dict[str, tuple[int, ...]] | None = None) -> dict:
    """Every check of a repacked file against its original. Returns problems (must be empty),
    expected differences and the storage of the compressed layers."""
    problems: list[str] = []
    expected: list[str] = []
    layers: dict[str, list] = {}
    if tg is None:
        tg = targets(orig)
    with h5py.File(orig, "r") as fa, h5py.File(new, "r") as fb:
        sa, sb = (
            fa.id.get_create_plist().get_version()[0],
            fb.id.get_create_plist().get_version()[0],
        )
        if sa != sb:
            problems.append(f"superblock version {sa} -> {sb}")
        if fa.userblock_size != fb.userblock_size:
            problems.append("userblock differs")
        names_a: list[str] = []
        names_b: list[str] = []
        fa.visit(names_a.append)
        fb.visit(names_b.append)
        if names_a != names_b:
            problems.append("object lists differ")
        for n in [""] + names_a:
            if n and n not in fb:
                continue
            a, b = (fa, fb) if not n else (fa[n], fb[n])
            if type(a) is not type(b):
                problems.append(f"/{n}: object type differs")
                continue
            cmp = attrs_compare(fa, a, fb, b)
            if cmp == "differ":
                problems.append(f"/{n}: attributes differ")
            elif cmp == "order":
                problems.append(f"/{n}: attribute order differs")
            if isinstance(a, h5py.Group):
                pa = h5py.h5g.open(fa.id, ("/" + n).encode()).get_create_plist()
                pb = h5py.h5g.open(fb.id, ("/" + n).encode()).get_create_plist()
                if pa.get_link_creation_order() != pb.get_link_creation_order():
                    problems.append(f"/{n}: link creation order flags differ")
                if pa.get_attr_creation_order() != pb.get_attr_creation_order():
                    problems.append(f"/{n}: group attribute creation-order flags differ")
                continue
            if isinstance(a, h5py.Datatype):
                if a.dtype != b.dtype:
                    problems.append(f"/{n}: committed datatype differs")
                continue
            if (a.shape, a.dtype, a.maxshape) != (
                b.shape,
                b.dtype,
                b.maxshape,
            ) or a.id.get_type() != b.id.get_type():
                problems.append(f"/{n}: shape, dtype or maxshape differs")
                continue
            if np.asarray(a.fillvalue).tobytes() != np.asarray(b.fillvalue).tobytes():
                problems.append(f"/{n}: HDF5 fill value differs")
            if (
                a.id.get_create_plist().get_attr_creation_order()
                != b.id.get_create_plist().get_attr_creation_order()
            ):
                problems.append(f"/{n}: dataset attribute creation-order flags differ")
            if h5py.h5ds.is_scale(a.id) != h5py.h5ds.is_scale(b.id):
                problems.append(f"/{n}: dimension-scale flag differs")
            if [[s.name for s in d.values()] for d in a.dims] != [
                [s.name for s in d.values()] for d in b.dims
            ] or [d.label for d in a.dims] != [d.label for d in b.dims]:
                problems.append(f"/{n}: dimension scales differ")
            if n in tg:
                want = (tg[n], "gzip", 4, True, False, None)
                if storage(a) != (None, None, None, False, False, None) or storage(b) != want:
                    problems.append(f"/{n}: storage {storage(a)} -> {storage(b)}, wanted {want}")
                layers[n] = [a.id.get_storage_size(), b.id.get_storage_size()]
            elif storage(a) != storage(b):
                problems.append(f"/{n}: storage changed {storage(a)} -> {storage(b)}")
            if not same_values(a, b):
                problems.append(f"/{n}: values differ")
        missing = sorted(set(tg) - set(layers))
        if missing:
            problems.append(f"layers not compressed: {missing}")

    # The netCDF view: the header is identical, and the special attributes (-s) differ only by the
    # storage of the compressed layers
    if ncdump(orig, "-h") != ncdump(new, "-h"):
        problems.append("ncdump -h differs")
    removed, added = [], []
    for line in difflib.unified_diff(ncdump(orig, "-hs"), ncdump(new, "-hs"), lineterm="", n=0):
        if line.startswith(("---", "+++", "@@")):
            continue
        (removed if line[0] == "-" else added).append(line[1:].strip())
    short = {n.rsplit("/", 1)[-1] for n in tg}
    for line in removed:
        var, _, rest = line.partition(":")
        if var not in short or rest != '_Storage = "contiguous" ;':
            problems.append(f"ncdump -hs: unexpected removed line {line!r}")
    for line in added:
        var, _, rest = line.partition(":")
        if var not in short or not NC_STORAGE_ATTRS.match(rest):
            problems.append(f"ncdump -hs: unexpected added line {line!r}")
    if len(removed) != len(tg) or len(added) != 4 * len(tg):
        problems.append(
            f"ncdump -hs: {len(removed)} lines removed and {len(added)} added for {len(tg)} layers"
        )
    return {"problems": problems, "expected_differences": sorted(set(expected)), "layers": layers}


def drop_cache(path: Path, sync: bool = False) -> None:
    fd = os.open(path, os.O_RDONLY)
    try:
        if sync:
            os.fsync(fd)
        os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
    finally:
        os.close(fd)


def sha256_direct(path: Path) -> str:
    """SHA-256 of a file read with direct I/O (from the disk or NFS server, not the page cache)."""
    dd = subprocess.Popen(
        ["dd", f"if={path}", "iflag=direct", "bs=16M", "status=none"], stdout=subprocess.PIPE
    )
    h = hashlib.sha256()
    for block in iter(lambda: dd.stdout.read(1 << 24), b""):
        h.update(block)
    if dd.wait() != 0:
        raise OSError(f"direct read failed: {path}")
    return h.hexdigest()


def keep_metadata(src: Path, dst: Path) -> None:
    st = src.stat()
    os.chmod(dst, st.st_mode & 0o7777)
    os.utime(dst, ns=(st.st_atime_ns, st.st_mtime_ns))


def place(src: Path, dst: Path) -> dict[str, Any]:
    """Repack (.nc) or copy (other files) src to dst, check it, and only then put it in place.
    Never overwrites dst."""
    part = dst.parent / f".{dst.name}.part"
    st = src.stat()
    rec: dict[str, Any] = {"src": str(src), "dst": str(dst), "src_bytes": st.st_size}
    t0 = time.time()
    try:
        if dst.exists():
            raise FileExistsError(f"{dst} already exists")
        dst.parent.mkdir(parents=True, exist_ok=True)
        part.unlink(missing_ok=True)
        if src.suffix == ".nc":
            rec["kind"] = "repack"
            tg = targets(src)
            rebuild(src, part, tg)
            drop_cache(part, sync=True)
            rec.update(check_pair(src, part, tg))
            drop_cache(part)
        else:
            rec["kind"] = "copy"
            subprocess.run(
                [
                    "dd",
                    f"if={src}",
                    f"of={part}",
                    "iflag=direct",
                    "oflag=direct",
                    "bs=16M",
                    "conv=fsync",
                    "status=none",
                ],
                check=True,
            )
            a, b = sha256_direct(src), sha256_direct(part)
            rec.update(sha256=b, problems=[] if a == b else [f"SHA-256 differs: {a} vs {b}"])
        if rec["problems"]:
            part.unlink(missing_ok=True)
            rec["status"] = "failed"
        else:
            keep_metadata(src, part)
            os.replace(part, dst)
            rec.update(status="ok", dst_bytes=dst.stat().st_size)
    except Exception as exc:
        part.unlink(missing_ok=True)
        rec.update(status="failed", problems=[f"{type(exc).__name__}: {exc}"])
    rec["seconds"] = round(time.time() - t0, 1)
    return rec


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__)
        return 2
    rec = place(Path(argv[0]), Path(argv[1]))
    print(json.dumps(rec))
    return 0 if rec["status"] == "ok" else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
