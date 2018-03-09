import numpy as np
import skimage.io
import paramiko
import frontmatter
import pandas as pd
import scp
import yaml
import os


def scrape_frontmatter(dirname, file='README.md'):
    """
    Reads the status of a given experimental dataset. This status is embedded
    in the README.md file as a YAML metadata block.

    Parameters
    ----------
    dirname : str
        Directory from which to parse.
    file: str
        Name of file containing YAML frontmatter. Default is 'README.md'

    Returns
    -------
    info : dict or pandas DataFrame
        A dictionary with all frontmatter keys and values.

    Raises
    ------
    UserWarning
        A UserWarning is raised if the scraped yaml frontmatter does not have
        a 'status' key or the value is not in `['accepted', 'rejected',
        'questionable']`.
    """
    # Grab file from directory.
    if dirname[-1] == '/':
        filename = '{}{}'.format(dirname, file)
    else:
        filename = '{}/{}'.format(dirname, file)

    # Scrape and return as desired.
    with open(filename) as f:
        info, _ = frontmatter.parse(f.read())
    if 'status' not in info.keys():
        raise UserWarning(
            'key `status` not found in metadata keys. Skipping {}'.format(dirname))

        info = {}
    elif info['status'].lower() not in ['accepted', 'questionable', 'rejected']:
        raise UserWarning('Value `status: {}` not an acceptable flag. Skipping {}'.format(
            info['status'].lower(), dirname))
        info = {}
    return info


def _create_SSHclient(username, port, server='delbruck'):
    """
    Establishes an SSH connection to one of the servers.

    Parameters:
    -----------
    username : str
        Username on the server
    port : int
        Super secret port number
    server : str
        First name of the beloved server.

    Returns
    -------
    client : SSH Object
        Client object connected to remote server.

    Notes
    -----
    This function is not compatabile with password-only login. You must have
    the proper system keys installed.

    Acknowledgements
    ----------------
    This function taken from a StackOverflow anwer by Tom Shen.
    https://stackoverflow.com/questions/250283/how-to-scp-in-python
    """
    client = paramiko.SSHClient()
    client.load_system_host_keys()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect('{}@caltech.edu'.format(server), port, username)
    return client


def yank_clist(src, username, port, dest, server='delbruck'):
    """
    Pulls the clist files from a target directory on delbrück to local.

    Parameters
    ----------
    src: str
        Relative path to the root diriectory for the experiment.
    dst: str
        Target directory on local. If None, the clists will be pulled into
        `data/images` with the same folder name.
    username: str
        Username to access the data.
    port : int
        Super secret port number.
    """
    prefix = 'git/mwc_growth/data/images/'

    # Avoid very stupid errors with the provided source directory.
    if src[-1] == '/':
        src = src[:-1]

    # Define the source pattern.
    src = '{}{}xy\*/clist.mat'.format(prefix, src)
    # Connect to the remote server.
    client = _create_SSHclient(username, port, server=server)
    scp_obj = scp.SCPClient(client.get_transport())

    # Transfer the files.
    scp_obj.get(src, dest)

    # Kill the connection
    client.close()
    print("clist files successfully transferred.")
