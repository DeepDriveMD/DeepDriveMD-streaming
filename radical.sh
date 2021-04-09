module load gcc/7.3.1

export PATH=$HOME/spack/bin:$PATH

. /etc/profile.d/conda.sh
conda activate /usr/workspace/cv_ddmd/conda/rct.devel

shopt -s direxpand

export IBM_POWERAI_LICENSE_ACCEPT=yes
export WORKSPACE=/usr/workspace/$USER
export SHARED=/usr/workspace/cv_ddmd/$USER
export GPFS1=/p/gpfs1/$USER
export SB=$GPFS1/radical.pilot.sandbox
export MyDir=$1
export MyDir1=$GPFS1/DDMD/$MyDir/entk_cvae_md
export MyDir2=$MyDir1/analysis
export PATH=$HOME/bin:$HOME/spack/bin:$PATH

shopt -s direxpand

alias wks='cd $WORKSPACE'
alias gpfs='cd $GPFS1'

alias sb='cd $SB'
alias sb1='cd $SB/$(ls -d $MyDir1/re.* | xargs basename)/pilot.0000'

alias cd1='cd $MyDir1'
alias cd2='cd $MyDir2'
alias cd3='cd $SHARED'

alias bu='bkill -u $USER'
alias bj='bjobs'
alias topU='top -U $USER'

function mcd1()
{
    cd $GPFS1/DDMD/$1/entk_cvae_md
}

function mcd2()
{
    cd $GPFS1/DDMD/$1/entk_cvae_md/analysis
}

function msb1()
{
    cd $SB/$(ls -d $GPFS1/DDMD/$1/entk_cvae_md/re.* | xargs basename)/pilot.0000
}

function em
{
    emacs -nw +$2:1 $1
}
