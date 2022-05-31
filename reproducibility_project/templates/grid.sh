{% extends "slurm.sh" %}

{% block header %}
{% set gpus = operations|map(attribute='directives.ngpu')|sum %}
    {{- super () -}}

{% set cpus = operations|map(attribute='directives.np')|sum %}
    {{- super () -}}

{% if gpus %}
#SBATCH -q gpu
#SBATCH --gres gpu:{{ gpus }}
{%- else %}
#SBATCH -q primary
#SBATCH --constraint=intel
{%- endif %}

#SBATCH -N 1
#SBATCH --ntasks {{ cpus }}
#SBATCH -o output-%j.dat
#SBATCH -e error-%j.dat

echo  "Running on host" hostname
echo  "Time is" date
source /sw/packages/python-tools/anaconda3/etc/profile.d/conda.sh
conda activate mamba


{% endblock header %}

{% block body %}
    {{- super () -}}


{% endblock body %}
