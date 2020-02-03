import os, subprocess, time


def run_command(command):
    myrun = subprocess.Popen(command, shell=True)
    pid = myrun.pid
    #print("Running {} as {} ".format(command, myrun))
    observe_job(process=myrun, justwait=True)
    return pid


def observe_job(process, justwait):
    if justwait:
        # process.wait()
        os.waitpid(process.pid, 0)
    else:
        #print("Actively monitoring ", process, process.poll(), process.pid)
        while check_pid(process.pid):
            time.sleep(1)
    #print("process {} is done".format(process))


def check_pid(pid):
    """ Check For the existence of a unix pid. Copied from stackoverflow."""
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    else:
        return True


def cleanup_dir(dirpath):
    files = os.listdir(dirpath)
    # add wildcard filters if needed
    for file in files:
        os.remove(os.path.join(dirpath, file))


def read_file(fname):
    with open(fname, "r") as fout:
        text = fout.read()
    return text


def write_file(fname, text):
    with open(fname, 'w') as fin:
        fin.write(text)


def prepare_dir(dirpath):
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)
