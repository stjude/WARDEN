#!/usr/bin/env python

import yaml
import sys
import subprocess
parameter_file = sys.argv[1]
# parameter_file="../test.yaml"
stream = file(parameter_file, 'r')


yaml_obj = (yaml.load(stream))


calc_norm_factors_method = yaml_obj["calc_norm_factors_method"]
contrasts = yaml_obj["contrasts"]
contrast_list = []
for c in contrasts:
    one, two = c.split("-")
    contrast_name = one + "_vs_" + two + "=" + c
    contrast_list.append(contrast_name)
contrast_string = ",".join(contrast_list)
filter_count = str(yaml_obj["filter_count"])
p_value_adjust = yaml_obj["p_value_adjust"]
filter_count = str(yaml_obj["filter_count"])
filter_count_type = yaml_obj["filter_count_type"]
# RPARAM=open("/R/R_parameters.txt",'w')
print("filter_count\t" + filter_count)
print("filter_count_type\t" + filter_count_type)
print("contrasts\t" + contrast_string)
print("calc_norm_factors_method\t" + calc_norm_factors_method)
print("p_value_adjust\t" + p_value_adjust)
# RPARAM.close()


if __name__ == '__main__':
    pass