#!/usr/bin/env python3
import asyncio
from itertools import chain

import utils
import sys
import subprocess
from datetime import timedelta
from pathlib import Path
from yapapi.log import enable_default_logger, log_summary, log_event_repr
from yapapi.runner import Engine, Task, vm
from yapapi.runner.ctx import WorkContext

from melatonin_analysis import get_output_filename

LOCAL_DATASETS = Path(__file__).parent.resolve().joinpath("datasets")
LOCAL_RECEIVED = Path(__file__).parent.resolve().joinpath("received")
LOCAL_RECEIVED_MOGGED = LOCAL_RECEIVED.joinpath("mogged")
SPREADSHEET_GLOB = "**/*.xlsx"
DATASET_XLSX_GLOB_FMT2 = "**/{group}/{name}.xlsx"
DATASET_XLSX_GLOB_GROUP_FMT1 = "**/{group}/*.xlsx"
DATASET_XLSX_GLOB_NAME_FMT1 = "**/{name}.xlsx"

GOLEM_WORK = "/golem/work"
ANALYSIS_SCRIPT_NAME = "melatonin_analysis.py"
ANALYSIS_SCRIPT = f"{GOLEM_WORK}/{ANALYSIS_SCRIPT_NAME}"
ANIMATED_FILENAME = "animated.gif"

# whitelists of datasets to run, filtered by group & name (empty list = any)
# group_whitelist = []
group_whitelist = ["ppd"]
name_whitelist = []
# name_whitelist = ["Dver002ur", "Dver004ur"]
# name_whitelist = ["Dver002ur"]


def find_whitelisted_datasets():
    """
    :return: The datasets per the group/name whitelists.
    """
    datasets = iter(())
    if len(group_whitelist) == 0:
        if len(name_whitelist) == 0:
            # all groups, all names
            datasets = chain(datasets, LOCAL_DATASETS.glob(SPREADSHEET_GLOB))
        else:
            # any group, specific names
            for n in name_whitelist:
                datasets = chain(datasets, LOCAL_DATASETS.glob(DATASET_XLSX_GLOB_NAME_FMT1.format(name=n)))
    else:
        if len(name_whitelist) == 0:
            # specific groups, any names
            for g in group_whitelist:
                datasets = chain(datasets, LOCAL_DATASETS.glob(DATASET_XLSX_GLOB_GROUP_FMT1.format(group=g)))
        else:
            # specific groups, specific names
            for g in group_whitelist:
                for n in name_whitelist:
                    datasets = chain(datasets, LOCAL_DATASETS.glob(DATASET_XLSX_GLOB_FMT2.format(group=g, name=n)))
    return datasets


async def main(args):
    package = await vm.repo(
        image_hash="3a7b022e165069bf24dec76f137681b492e24f51cde95ddcc8b286b1",
        min_mem_gib=0.5,
        min_storage_gib=2.0,
    )

    async def worker(ctx: WorkContext, tasks):
        async for task in tasks:
            dataset, threshold = task.data
            dataset_name = dataset.stem
            dataset_name_ext = dataset.name
            remote_dataset = f"{GOLEM_WORK}/{dataset_name_ext}"
            output_filename = get_output_filename(dataset_name, threshold)
            output_path_local = LOCAL_RECEIVED.joinpath(output_filename)
            ctx.send_file(ANALYSIS_SCRIPT_NAME, ANALYSIS_SCRIPT)
            ctx.send_file(str(dataset), remote_dataset)
            ctx.run("chmod", "a+x", ANALYSIS_SCRIPT)
            ctx.run("/bin/sh", "-c", f"{ANALYSIS_SCRIPT} {remote_dataset} {threshold}")
            ctx.download_file(f"{GOLEM_WORK}/{output_filename}", str(output_path_local))
            yield ctx.commit()
            task.accept_task(result=output_filename)

    async with Engine(
        package=package,
        max_workers=args.number_of_providers,
        budget=10.0,
        # timeout should be number of datasets / number of providers dependent
        timeout=timedelta(minutes=25),
        subnet_tag=args.subnet_tag,
        event_emitter=log_summary(log_event_repr),
    ) as engine:
        # setup several threshold analyses per dataset
        cmdline_params = []
        datasets = find_whitelisted_datasets()
        thresholds = (25, 50, 75)
        for dataset in datasets:
            for threshold in thresholds:
                cmdline_params.append((dataset, threshold))

        async for task in engine.map(worker, [Task(data=params) for params in cmdline_params]):
            print(
                f"{utils.TEXT_COLOR_CYAN}"
                f"Task computed: {task}, result: {task.output}"
                f"{utils.TEXT_COLOR_DEFAULT}"
            )


if __name__ == "__main__":
    parser = utils.build_parser("melatonin analysis")
    parser.add_argument("--number-of-providers", dest="number_of_providers", type=int, default=4)
    args = parser.parse_args()
    enable_default_logger(log_file=args.log_file)
    sys.stderr.write(
        f"Using subnet: {utils.TEXT_COLOR_YELLOW}{args.subnet_tag}{utils.TEXT_COLOR_DEFAULT}\n"
    )

    loop = asyncio.get_event_loop()
    task = loop.create_task(main(args))

    try:
        loop.run_until_complete(task)

        # create animated .gif
        print("Mogrifying images...")
        subprocess.run(["mogrify", "-path", str(LOCAL_RECEIVED_MOGGED), "-format", "gif", "-resize", "50%", "*.png"],
                       cwd=LOCAL_RECEIVED)
        print("Creating animated .gif...")
        subprocess.run(["convert", "-delay", "20", "-loop", "0", "*.gif", ANIMATED_FILENAME],
                       cwd=LOCAL_RECEIVED_MOGGED)
        print(f"Complete.  Output: {LOCAL_RECEIVED_MOGGED.joinpath(ANIMATED_FILENAME)}")
    except (Exception, KeyboardInterrupt) as e:
        print(e)
        task.cancel()
        loop.run_until_complete(task)
