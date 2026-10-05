"""Command-line interface: argument parsing and routing to historic/forward processing."""
import argparse
import sys

from aria_coseis.modes import main_forward, main_historic


def main() -> None:
    """
    Run the main function based on the input arguments provided, in either 'historic' or 'forward' processing mode.
    Historic processing can be done for a single date or range of dates, with the option to specify the SLC pairing mode.
    Forward processing is used to generate co-seismic displacement products for new earthquakes.
    Example usage for historic processing: 
      python coseis_sar.py --historic --dates 2021-08-14 --pairing all
      python coseis_sar.py --historic --dates 2021-08-14 2021-09-07 --pairing all
      python coseis_sar.py --historic --dates 2014-06-14 2025-02-12 --pairing coseismic (all coseismic pairs from beginning of S1 data to 2025-02-12)
      python coseis_sar.py --historic --dates 2014-06-14 2025-02-12 --pairing coseismic --job_list --resolution 30 (only produce the job list for HYP3 processing, don't run any jobs locally) 
    Example usage for forward processing: 
      python coseis.py --forward
    """
    parser = argparse.ArgumentParser(
        description="Run historic, forward, or custom-list processing based on input arguments."
    )

    # Use a mutually exclusive group so users must pick only processing mode
    run_group = parser.add_mutually_exclusive_group(required=True)
    run_group.add_argument("--historic", action="store_true", help="Run historic processing.")
    run_group.add_argument("--forward", action="store_true", help="Run forward processing.")
    run_group.add_argument("--eq_list", type=str, help="Path to a custom JSON list of earthquakes to process.")

    parser.add_argument("--dates", nargs="+", help="Provide one or two dates in YYYY-MM-DD format for historic processing.")
    parser.add_argument("--aoi", help="Specify a path to a json file representing the area of interest (AOI).")
    parser.add_argument("--pairing", choices=["all", "sequential", "coseismic"], help="Specify the SLC pairing mode. Required for SAR processing.")
    parser.add_argument("--job_list", action="store_true", help="Create a list of jobs in HYP3 format for cloud processing.")
    parser.add_argument("--resolution", type=int, default=30, help="Output resolution for topsApp processing in meters. Default is 30m.")
    parser.add_argument("--sensor", choices=["sar", "sentinel-2", "landsat"], default="sar", 
                        help="Sensor: 'sar' (Sentinel-1), 'sentinel-2', or 'landsat'. Default is sar.")
    parser.add_argument("--do_processing", action="store_true", help="Execute local topsApp processing.")
    parser.add_argument("--send_email", action="store_true", help="Send email notifications.")
    parser.add_argument("--process_only", action="store_true", help="Skip discovery; only process existing jobs in the tracker.")
    parser.add_argument("--optical_backend", choices=["copernicus", "element84", "gee"], default="copernicus", help="Specify the optical data provider if sensor is optical. Default is copernicus.")
    parser.add_argument("--optical_level", choices=["raw", "toa", "sr"], default="toa", help="Specify the optical data level ('raw', 'toa', 'sr'). Default is 'toa'. Sentinel-2 does not support 'raw'.")

    args = parser.parse_args()

    # Global constraint check
    if args.sensor == 'sar' and ('--optical_backend' in sys.argv or '--optical_level' in sys.argv):
        print("Error: --optical_backend and --optical_level can only be used when --sensor is 'sentinel-2' or 'landsat'.")
        parser.print_help()
        exit(1)
        
    if args.sensor == 'sar' and not args.pairing:
        print("Error: --pairing is required when using --sensor sar. Options: 'all', 'sequential', 'coseismic'.")
        parser.print_help()
        exit(1)

    # Mode routing
    if args.historic:
        if args.send_email or args.do_processing:
            print("Error: --send_email and --do_processing are only supported in --forward mode.")
            exit(1)
        if not args.dates:
            print("Error: --dates is required when using --historic mode.")
            exit(1)
        if len(args.dates) > 2:
            print("Error: --dates should have at most two values (start_date [end_date]).")
            exit(1)

        start_date = args.dates[0]
        end_date = args.dates[1] if len(args.dates) == 2 else None
        main_historic(
            start_date=start_date, 
            end_date=end_date, 
            aoi=args.aoi, 
            pairing_mode=args.pairing, 
            job_list=args.job_list, 
            resolution=args.resolution, 
            sensor=args.sensor, 
            optical_backend=args.optical_backend,
            optical_level=args.optical_level
        )

    elif args.eq_list:
        if args.send_email or args.do_processing:
            print("Error: --send_email and --do_processing are only supported in --forward mode.")
            exit(1)
        if args.dates:
            print("Error: --dates cannot be used with --eq_list. The dates are defined in the file.")
            exit(1)
            
        main_historic(
            eq_list_path=args.eq_list, 
            aoi=args.aoi, 
            pairing_mode=args.pairing, 
            job_list=args.job_list, 
            resolution=args.resolution, 
            sensor=args.sensor, 
            optical_backend=args.optical_backend,
            optical_level=args.optical_level
        )

    elif args.forward:
        if args.sensor != 'sar':
            print("Error: --forward currently supports only --sensor sar.")
            exit(1)

        if args.job_list or args.dates or args.aoi:
            print("Error: --job_list, --dates, and --aoi cannot be used with --forward mode.")
            exit(1)
            
        if not args.do_processing and not args.send_email:
            print("Warning: Running --forward without --do_processing or --send_email. The script will only update tracking files.")

        if args.job_list:
            print("Error: --job_list is only supported in --historic mode.")
            parser.print_help()
            exit(1)
            
        if args.dates:
            print("Error: --dates cannot be used with --forward mode.")
            parser.print_help()
            exit(1)

        if not args.pairing:
            print("Error: --pairing is required when using --forward mode. Options: 'all', 'sequential', 'coseismic'.")
            parser.print_help()
            exit(1)

        if args.aoi:
            print("Error: --aoi cannot be used with --forward mode.")
            parser.print_help()
            exit(1)
            
        if not args.do_processing and not args.send_email:
            print("Warning: Running --forward without --do_processing or --send_email. The script will only update tracking files.")

        main_forward(
            pairing_mode=args.pairing, 
            resolution=args.resolution, 
            do_processing=args.do_processing, 
            send_email_flag=args.send_email, 
            process_only=args.process_only
        )
