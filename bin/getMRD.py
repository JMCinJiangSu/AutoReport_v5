#-*- coding:gbk -*-
import copy, re
from libs.rule import decimal_float

'''
Discription
	
	该脚本用来获取MRD结果。
'''

def getMRD(jsonDict):
	mrd_list = copy.deepcopy(jsonDict["mrd"]) if "mrd" in jsonDict.keys() and jsonDict["mrd"] else []
	for i in mrd_list:
		i["ctdna_conc"] = decimal_float(i["ctdna_conc"])
		i["receive_date"] = re.split(" ", i["receive_date"])[0] if i["receive_date"] else ""
	mrd_list = sorted(mrd_list, key = lambda i : i["receive_date"])

	mrd_last = mrd_list[-1] if mrd_list else []
	mrd = {
		"mrd_list" : mrd_list, 
		"mrd_last" : mrd_last
		}

	return mrd