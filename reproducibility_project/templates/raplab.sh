{% extends "slurm.sh" %}

{% block header %}
{% set gpus = operations|map(attribute='directives.ngpu')|sum %}
    {{- super () -}}

{% if gpus %}
#SBATCH --partition gpu
#SBATCH --gres gpu:{{ gpus }}
{%- else %}
#SBATCH --partition gpu
{%- endif %}

#SBATCH -N 1
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
