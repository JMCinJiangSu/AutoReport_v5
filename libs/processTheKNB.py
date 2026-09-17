#-*- coding:gbk -*-
from libs.getEvi import varRegimen
import copy
from libs.rule import S_function

'''
Discription
	
	处理knb格式。 

'''

def process_knb(jsonDict, config):
	knb = copy.deepcopy(jsonDict["knb"])[0] if jsonDict["knb"] else {}
	# 2026.08.28-新增规则
	# 证据中有A/B级证据，KNB不显示
	if knb:
		evi_level_list = [evi["evi_conclusion"][0] for evi in knb["evi_sum"]]
		if set(["A", "B"]) & set(evi_level_list):
			knb["evi_sum"] = varRegimen(jsonDict, knb["evi_sum"], config, knb)
			clinic_num_s, knb["top_level"] = S_function(knb)
		else:
			return {}

	return knb