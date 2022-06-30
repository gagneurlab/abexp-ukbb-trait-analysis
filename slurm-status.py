#!/usr/bin/env python3
import re
import subprocess as sp
import shlex
import sys
import time
import logging
import argparse

logger = logging.getLogger("__name__")


parser = argparse.ArgumentParser(description='Get the status of some slurm job.')
parser.add_argument(
    '--status_attempts',
    action='store',
    dest="status_attempts",
    type=int,
    default=20,
    help='number of attenpts to obtain cluster status'
)
parser.add_argument(
    '--cluster',
    action='store',
    dest="cluster",
    default="",
    help='cluster name'
)

parser.add_argument(
    '--debug',
    action='store_true',
    dest="debug",
    help='Change log level to debug'
)

parser.add_argument(
    'jobid',
    action='store',
    help='Slurm job id'
)

args = parser.parse_args()

if args.debug is True:
    loglevel = logging.DEBUG
else:
    loglevel = logging.WARNING

jobid = args.jobid
STATUS_ATTEMPTS = args.status_attempts
if args.cluster != "":
    cluster = f"--cluster={args.cluster}"
else:
    cluster = ""

# setup logging
logging.basicConfig(
    level=loglevel,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("slurm-status.log"),
        # logging.StreamHandler()
    ]
)


def run_scontrol(cluster, jobid):
    sctrl_res = sp.check_output(
        shlex.split(f"scontrol {cluster} -o show job {jobid}")
    )
    m = re.search(r"JobState=(\w+)", sctrl_res.decode())
    res = {jobid: m.group(1)}

    if res[jobid] == "PREEMPTED":
        m = re.search(r"Requeue=(\w+)", sctrl_res.decode())
        requeueable = m.group(1)

        if requeueable == "1":
            res[jobid] = "REQUEUED"

    return res


def run_sacct(cluster, jobid):
    sacct_res = sp.check_output(shlex.split("sacct -P -b -j {} -n".format(jobid)))
    res = {x.split("|")[0]: x.split("|")[1] for x in sacct_res.decode().strip().split("\n")}
    return res


for i in range(STATUS_ATTEMPTS):
    try:
        res = run_sacct(cluster, jobid)
        break
    except sp.CalledProcessError as e:
        logger.error("sacct process error")
        logger.error(e)
    except IndexError as e:
        pass

    # Try getting job with scontrol instead in case sacct is misconfigured
    try:
        res = run_scontrol(cluster, jobid)
        break
    except sp.CalledProcessError as e:
        logger.error("scontrol process error")
        logger.error(e)
        if i >= STATUS_ATTEMPTS - 1:
            print("failed")
            exit(0)
        else:
            time.sleep(1)

status = res[jobid]
logger.info("%s: '%s'", jobid, status)

if status == "BOOT_FAIL":
    print("failed")
elif status == "OUT_OF_MEMORY":
    print("failed")
elif status.startswith("CANCELLED"):
    print("failed")
elif status == "COMPLETED":
    print("success")
elif status == "DEADLINE":
    print("failed")
elif status == "FAILED":
    print("failed")
elif status == "NODE_FAIL":
    print("failed")
elif status == "REQUEUED":
    print("running")
elif status == "PREEMPTED":
    print("failed")
elif status == "TIMEOUT":
    print("failed")
# Unclear whether SUSPENDED should be treated as running or failed
elif status == "SUSPENDED":
#    print("failed")
    print("running")
else:
    print("running")
