#-*- coding:gbk -*-
import copy,re
from libs import listResultToDict
from libs.getEvi import varRegimen
from libs.rule import S_function

'''
Discription
	
	该脚本用来获取MSI结果。
	提取内容：
	1. 结果
	2. 图片
	3. 治疗方案
	用于MSI特殊处理内容 

'''
# 识别是否为数值
def is_number(i):
	try:
		float(i)
		return True
	except:
		pass 
	if i.isnumeric():
		return True
	return False

def getMSI(jsonDict, config):
	msi_dict = listResultToDict.ListToDict(copy.deepcopy(jsonDict["msi"]))
	# MSI_num json无返回，通过脚本计算
	if msi_dict:
		# 2025.04.16-msi_score兼容百分数
		msi_dict["msi_score"] = float(msi_dict["msi_score"].replace("%", ""))/100 if re.search("%", str(msi_dict["msi_score"])) else float(msi_dict["msi_score"]) if msi_dict["msi_score"] and is_number(msi_dict["msi_score"]) else 0
		# 2025.04.16-兼容完成
		msi_dict["msi_num_cp40"] = int(float(msi_dict["msi_score"]) * 55 + 0.5) if msi_dict["msi_score"] or msi_dict["msi_score"]==0 else ""
		# 2024.12.13-新增一个116的MSI_num，并且将msi_socre转化为百分数
		msi_dict["msi_num_116"] = int(float(msi_dict["msi_score"]) * 53 + 0.5) if msi_dict["msi_score"] or msi_dict["msi_score"]==0 else ""
		msi_dict["msi_score_str"] = "{:.2%}".format(float(int(float(msi_dict["msi_score"]) * 10000 + 0.5) / 10000)) if msi_dict["msi_score"] else 0
		# 2024.12.13-新增完成

	if "evi_sum" in msi_dict.keys() and msi_dict["evi_sum"]:
		msi_dict["evi_sum"] = varRegimen(jsonDict, msi_dict["evi_sum"], config, msi_dict)
		clinic_num_s, msi_dict["top_level"] = S_function(msi_dict)
	return msi_dict