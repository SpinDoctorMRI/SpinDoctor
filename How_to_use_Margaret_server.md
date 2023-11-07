# How to use Margaret server

The server is connected via ssh service. 
For Windows 10 and later, the ssh service is built in.
You can generate new key by the command `ssh-keygen`.

Then you can find two new files:

 - ~/.ssh/id_rsa.pub
  - ~/.ssh/id_rsa (private key)

P.S. Never give the private key to the others.

## Ask the permission

To get access to Margaret server, you need to open a ticket to ask the permission:

1. Login to the Inria's intranet, click `helpdesk`
2. Click `Soumettre une demande informatique`, then `Demande de service`, then `Demande sur moyens de calcul`
3. Write the demande. It requires (1) the machine name/computer number, and also (2) the ~/.ssh/id_rsa.pub. 

Once the staff adds you computer into the server trusted list, you can connect to the server.

## First time configuration

1. Create a new ssh configuration file (if it does not exist) `~/.ssh/config`, and configure with

    ~~~
    # contents of ~/.ssh/config
    Host margaret
    HostName margaret.saclay.inria.fr
    User YourUserName
    ProxyJump YourUserName@ssh.saclay.inria.fr

    Host marg*
    HostName %h
    ProxyJump margaret
    User YourUserName
    ~~~

2. Open a terminal, tap `ssh margaret`. You should normally login into the portal node of the server. You should see

    ~~~
    [YourUserName@margaret ~]$
    ~~~

    in terminal.


    The Margaret server is organized by Slurm. The node 'margaret' is only a jump node. You can manage and transfer files or submit job under this node. But do not run the computation directly on this node. You should run the computation on compute nodes, for example 'marg008', 'margpu010', etc.

3. In order to connect the compute nodes directly from your computer, You need to add the local ssh key to the server, by the command

    `cat ~/.ssh/id_rsa.pub | ssh margaret 'cat - >> ~/.ssh/authorized_keys'`

4. Create directory under `/scratch` by the command 

    `mkdir /scratch/$USER`

    On the server, the directory `/home/idefix/YourUserName` is only for config files. The `/home/` section is quite small, so you should put the codes and data, and install Python under `/scratch/YourUserName/`.

    P.S. If you want bigger space, you need to open a new ticket and ask for the access to the `/data3` section

5. Install conda environment by the command 

    ~~~
    cd /scratch/YourUserName
    BASE_URL=https://github.com/conda-forge/miniforge/releases/latest/download/
    wget $BASE_URL/Mambaforge-Linux-x86_64.sh

    bash Mambaforge-Linux-x86_64.sh -b -p /scratch/YourUserName/mambaforge

    # Activate conda env by default
    echo 'source /scratch/YourUserName/mambaforge/etc/profile.d/conda.sh
    conda activate
    ' >> $HOME/.bashrc
    source $HOME/.bashrc
    ~~~

6. Create new environment with lastest Python version by the command

    `conda create -n NameOfEnv python=3.10`


7. Require packages for running the codes for simulation-driven framework paper:

Numpy, matplotlib, pandas, nibabel, nilearn, seaborn, jupyter, plotly, scipy, scikit-learn, dipy, python-spams, tqdm, joblib, setuptools: `conda install PackageName`

If one package does not exist in conda, you can try pip. For example, `pip install spams`

Amico: `pip install dmri-amico`

Pytorch: `conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch-nightly -c nvidia`


## Regular routine to connect to Margaret using VS Code

1. Open a terminal, tap `ssh margaret`. Once connected, tap one of the following commands:

    - `srun --pty -c 16 bash` to connect to a CPU node with 16 cores, with interactive session, named bash;

    - `srun --pty -c 128 bash` to connect to a CPU node with 128 cores, with interactive session, named bash;

    - `ssh margaret -t srun --pty --partition=gpu --gres=gpu:1 bash` to connect to a GPU node, with interactive session, named bash;

    - `srun --pty -c 16 --time 4:00:00 bash` to connect to a CPU node with 16 cores, with interactive session, named bash. And it will be closed after 4 hours. This one is useful because during maintenance, all the jobs should have an ending time.

    Then Slurm will assign one compute node, for example, 'marg008'.

    Keep this terminal alive.

2. Open VS Code, tap `F1`, then `Remote-SSH: Connect to host`, then tap `marg008`. Normally you should connect to marg008 successfully. You can see the node's name on the left-bottom corner.

3. Choose `File`, then `Open folder`, then tap `/scratch/YourUserName`. The working directory is changed to `/scratch/YourUserName`. You can easily open and edit files of this directory. And this is an interactive session, you can run the computation directly (like Python scripts, jupyter notebooks etc.)

4. To run python scripts on CPU node, create a `job.sh` file with
    ~~~
    #!/bin/bash

    #SBATCH --job-name=idefix
    #SBATCH --output=output_%j.txt
    #SBATCH --error=error_%j.txt

    #SBATCH --time=24:00:00
    #SBATCH --ntasks=64 ## number of nodes required

    ## load modules
    source /scratch/YourUserName/mambaforge/bin/activate NameOfEnv

    ## execution
    python path_of_the_file.py
    ~~~
    and open a new terminal in VS Code. Change the directory of this terminal to the location of the `job.sh` script. Run the command

    `sbatch job.sh`

5. To run python scripts on GPU node, it is similar but change the script to
    ~~~
    #!/bin/bash

    #SBATCH --job-name=idefix
    #SBATCH --output=output_%j.txt
    #SBATCH --error=error_%j.txt

    #SBATCH --time=24:00:00
    #SBATCH --partition=gpu
    #SBATCH --gres=gpu:1

    ## load modules
    source /scratch/YourUserName/mambaforge/bin/activate NameOfEnv

    ## execution
    python path_of_the_file.py
    ~~~

6. To run the Matlab scripts on CPU node, change the script to
    ~~~
    #!/bin/bash

    #SBATCH --job-name=idefix
    #SBATCH --output=output_%j.txt
    #SBATCH --error=error_%j.txt

    #SBATCH --time=24:00:00
    #SBATCH --ntasks=64 ## number of nodes required

    ## load modules
    module load matlab

    ## execution
    matlab -batch "path_of_the_file"
    ~~~
    


## Some useful Slurm commands

 - `sinfo`, introspect the state of nodes
 - `squeue -u YourUserName`, introspect the state of your jobs
 - `scancel JobID`, cancel your job


## Useful link
The link about using the margaret server:
https://docs.google.com/presentation/d/1eLhb5vbrM2menXzrAxqzeZ_p9msNWxjWtGiasYZHaf4/edit#slide=id.gef73fd5f03_0_284

# Transfer files with server

Personally I use WinSCP on Windows and on Mac.


Configuration is quite simple. The File protocol is SFTP, Host name is `margaret.saclay.inria.fr`, Port number is 22, User name is the Inria account name.