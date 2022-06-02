{% extends "slurm.sh" %}

{% block header %}
{{- super () -}}
{% set gpus = operations|map(attribute='directives.ngpu')|sum %}
{% set cpus = operations|map(attribute='directives.np')|sum %}

{% if gpus %}
#SBATCH -q gpu
#SBATCH --gres gpu:{{ gpus }}
{%- else %}
#SBATCH -q primary
{%- endif %}

#SBATCH -N 1
#SBATCH --ntasks {{ cpus }}
#SBATCH -o output-%j.dat
#SBATCH -e error-%j.dat

echo  "Running on host" hostname
echo  "Time is" date
source ~/anaconda3/etc/profile.d/conda.sh
conda activate mosdef-study38
module load gnu
{% endblock header %}

{% block body %}
    {{- super () -}}


{% endblock body %}
