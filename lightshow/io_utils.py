from subprocess import Popen, PIPE
import time


def run_command(cmd):
    """Execute the external command and get its exitcode, stdout and
    stderr."""

    t0 = time.time()
    proc = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
    out, err = proc.communicate()
    exitcode = proc.returncode
    dt = time.time() - t0

    return {
        "exitcode": exitcode,
        "output": out.decode("utf-8").strip(),
        "error": err.decode("utf-8").strip(),
        "elapsed": dt,
    }
