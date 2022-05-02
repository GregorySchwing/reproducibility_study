{# Templated in accordance with: https://www.olcf.ornl.gov/for-users/system-user-guides/summit/running-jobs/ #}
{% extends "lsf.sh" %}
{% block tasks %}
    {% set threshold = 0 if force else 0.9 %}
    {% set nn = operations|map('guess_resource_sets')|calc_num_nodes(parallel) %}
#BSUB -nnodes {{ nn }}
{% endblock tasks %}
{% block header %}
    {{- super() -}}
    {% set account = account|default(environment|get_account_name, true) %}
    {% if account %}
#BSUB -P {{ account }}
    {% endif %}

echo  "Running on host" hostname
echo  "Time is" date
source /sw/packages/python-tools/anaconda3/etc/profile.d/conda.sh
conda activate mamba

{% endblock header %}

{% block body %}
{{- super () -}}


{% endblock body %}
