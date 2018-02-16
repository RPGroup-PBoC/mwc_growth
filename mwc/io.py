import numpy as np
import skimage.io
import paramiko
import scp
import yaml
import os


def fetch_status(dir, file='README.md'):
    """
    Reads the status of a given experimental dataset. This status is embedded
    in the README.md file as a YAML metadata block.
    """
    return True


def pull_clist(src, username, port, dest):
    """
    Pulls the clist files from a target directory on delbrück to local.

    Parameters
    ----------
    src: str
        Source directory on Delbrück. This should be given as a relative path.
    dst: str
        Target directory on local. If None, the clists will be pulled into
        `data/images` with the same folder name.
    username: str
        Username to access the data.
    port : int
        Super secret Delbrück port number.
    """
    return True


def deoposit_data(src, dest, port, username):
    return True
