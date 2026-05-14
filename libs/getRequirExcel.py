#-*- coding:gbk -*-
import os
import re
import xlrd

'''
Discription
	
	获取config文件夹中的报告模板配置文件，转化为便于其他程序调用的格式。 
	{
		临检 : {
			定制 : {
				*新增 (医院，产品1，项目) ：模板名1,
				(医院, 产品1) : 模板名1
				},
			通用-简版 : {
				产品1 : 模板名1
				},
			通用-完整 : {
				产品1 : 模板名1
				}
			},
		进院 : {
			定制 : {
				(医院, 科室, 产品1) : 模板名1,
				(医院, 产品1) : 模板名1
 				},
			通用 : {
				产品1 : 模板名1
				}
		}
	
'''

#BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#requirenment_path = os.path.join(BASE_DIR, "config/report_requirenment.xlsx")
# 2024.12.31-新增CNV类型，用于模板判断展示MLPA还是CNV

key_stran = {
	"产品名称" : "prod_name",
	"检测业务类型" : "business_type",
	"模板类型" : "report_type",
	"单位名称" : "company",
	"科室" : "hosp_depart",
	"模板名" : "report_name",
	"状态" : "status",
	"模板开发者" : "developer",
	"添加人" : "auditors",
	"添加时间" : "auditors_time",
	"备注" : "note",
	"更新记录" : "update",
	"临检" : "rummage",
	"进院" : "hospital",
	"定制" : "CustomEdition",
	"通用" : "Universal",
	"通用-简版" : "Universal_simple",
	"通用-完整" : "Universal_complete",
	"基础模板拼接" : "merge_template",
	"基础模板-动态" : "temp_dynamic",
	"基础模板-固定" : "temp_fixed",
	"药企项目名称" : "product_name",
	"CNV类型" : "judge_brca_cnv",
	"订单类型" : "order_type",
	"免费类型" : "free_type",
	"来源" : "order_origin"
	}

def getRequirenment(config):
	#BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
	requirenment_path = os.path.join(config, "report_requirenment.xlsx")
	xls = xlrd.open_workbook(requirenment_path)
	requirenment_sheet = xls.sheet_by_name("report_requirenment")
	key = requirenment_sheet.row_values(0)
	Data = []
	for num in range(1, requirenment_sheet.nrows):
		rows = requirenment_sheet.row_values(num)
		tmpdict = {}
		for i in range(len(key)):
			tmpdict[key_stran.get(key[i])] = key_stran[rows[i]] if rows[i] in key_stran.keys() else rows[i]
		Data.append(tmpdict)
	
	# 2025.09.29-兼容多产品出一份报告，配置中产品名称填写多个，按，分隔，这边需要重新排序下
	for i in Data:
		i["prod_name"] = "，".join(sorted(re.split("，", i["prod_name"])))
	# 2025.09.29-兼容完成

	#配置信息字典结构	
	requir_dict = {
		"rummage" : {
			"CustomEdition" : {},
			"Universal_simple" : {},
			"Universal_complete" : {}
			},
		"hospital" : {
			"CustomEdition" : {},
			"Universal" : {}
			}
		}
	
	# 新增BRCA CNV字典结构-2024.12.31
	brca_cnv_dict = {
		"rummage" : {
			"CustomEdition" : {},
			"Universal_simple" : {},
			"Universal_complete" : {}
			},
		"hospital" : {
			"CustomEdition" : {},
			"Universal" : {}
			}
		}

	for i in Data:
		# 2026.01.21-新增来源限制-不为LIMS
		#if str(int(i["status"])) == "0":
		if str(int(i["status"])) == "0" and i["order_origin"] != "LIMS":
		# 2026.01.21-新增完成
			if i["company"]:
				if i["hosp_depart"]:
					requir_dict[i["business_type"]][i["report_type"]].setdefault((i["company"],i["hosp_depart"], i["prod_name"]), "")
					requir_dict[i["business_type"]][i["report_type"]][(i["company"], i["hosp_depart"], i["prod_name"])] = i["report_name"]
					brca_cnv_dict[i["business_type"]][i["report_type"]].setdefault((i["company"],i["hosp_depart"], i["prod_name"]), "")
					brca_cnv_dict[i["business_type"]][i["report_type"]][(i["company"], i["hosp_depart"], i["prod_name"])] = i["judge_brca_cnv"]

				else:
					# 新增药企配置-2023.11.17
					if i["product_name"]:
						requir_dict[i["business_type"]][i["report_type"]].setdefault((i["company"], i["prod_name"], i["product_name"]), "")
						requir_dict[i["business_type"]][i["report_type"]][(i["company"], i["prod_name"], i["product_name"])] = i["report_name"]
						brca_cnv_dict[i["business_type"]][i["report_type"]].setdefault((i["company"], i["prod_name"], i["product_name"]), "")
						brca_cnv_dict[i["business_type"]][i["report_type"]][(i["company"], i["prod_name"], i["product_name"])] = i["judge_brca_cnv"]
					else:
					# 新增结束-2023.11.17
						requir_dict[i["business_type"]][i["report_type"]].setdefault((i["company"], i["prod_name"]), "")
						requir_dict[i["business_type"]][i["report_type"]][(i["company"], i["prod_name"])] = i["report_name"]
						brca_cnv_dict[i["business_type"]][i["report_type"]].setdefault((i["company"], i["prod_name"]), "")
						brca_cnv_dict[i["business_type"]][i["report_type"]][(i["company"], i["prod_name"])] = i["judge_brca_cnv"]
			else:
					requir_dict[i["business_type"]][i["report_type"]].setdefault(i["prod_name"], "")
					requir_dict[i["business_type"]][i["report_type"]][i["prod_name"]] = i["report_name"]
					brca_cnv_dict[i["business_type"]][i["report_type"]].setdefault(i["prod_name"], "")
					brca_cnv_dict[i["business_type"]][i["report_type"]][i["prod_name"]] = i["judge_brca_cnv"]

	# 是否需要拼接基础模板对照表
	report_merge_dict = {}
	for i in Data:
		if str(int(i["status"])) == "0":
			if i["report_name"] not in report_merge_dict.keys():
				report_merge_dict.setdefault(i["report_name"],{})
			report_merge_dict[i["report_name"]] = {
				"merge_template" : str(int(i["merge_template"])),
				"temp_dynamic" : re.split(";", i["temp_dynamic"].replace("\n","").replace(" ","")) if i["temp_dynamic"] else "",
				"temp_fixed" : re.split(";", i["temp_fixed"].replace("\n","").replace(" ","")) if i["temp_fixed"] else ""
			}
	#print (requir_dict)

	return requir_dict, report_merge_dict, brca_cnv_dict

# 2026.01.21-新增v2-适配临检LIMS来源订单
def getRequirenment_v2(config):
	requirenment_path = os.path.join(config, "report_requirenment.xlsx")
	xls = xlrd.open_workbook(requirenment_path)
	requirenment_sheet = xls.sheet_by_name("report_requirenment")
	key = requirenment_sheet.row_values(0)
	Data = []
	for num in range(1, requirenment_sheet.nrows):
		rows = requirenment_sheet.row_values(num)
		tmpdict = {}
		for i in range(len(key)):
			tmpdict[key_stran.get(key[i])] = key_stran[rows[i]] if rows[i] in key_stran.keys() else rows[i]
		Data.append(tmpdict)
	
	# 2025.09.29-兼容多产品出一份报告，配置中产品名称填写多个，按，分隔，这边需要重新排序下
	for i in Data:
		i["prod_name"] = "，".join(sorted(re.split("，", i["prod_name"])))
	# 2025.09.29-兼容完成

	#配置信息字典结构	
	requir_dict = {
		"rummage" : {
			"CustomEdition" : {},
			"Universal_simple" : {},
			"Universal_complete" : {}
			},
		"hospital" : {
			"CustomEdition" : {},
			"Universal" : {}
			}
		}
	
	# 新增BRCA CNV字典结构-2024.12.31
	brca_cnv_dict = {
		"rummage" : {
			"CustomEdition" : {},
			"Universal_simple" : {},
			"Universal_complete" : {}
			},
		"hospital" : {
			"CustomEdition" : {},
			"Universal" : {}
			}
		}

	for i in Data:
		# 这边需要再加上来源不为BIMS
		if str(int(i["status"])) == "0" and i["order_origin"] != "BIMS":
			if i["company"]:
				# 先处理药企的
				if i["product_name"]:
					requir_dict[i["business_type"]][i["report_type"]].setdefault((i["company"], i["prod_name"], i["product_name"]), "")
					requir_dict[i["business_type"]][i["report_type"]][(i["company"], i["prod_name"], i["product_name"])] = i["report_name"]
					brca_cnv_dict[i["business_type"]][i["report_type"]].setdefault((i["company"], i["prod_name"], i["product_name"]), "")
					brca_cnv_dict[i["business_type"]][i["report_type"]][(i["company"], i["prod_name"], i["product_name"])] = i["judge_brca_cnv"]
				else:
					# 存在order_type（临检）
					if i["order_type"]:
						# 存在免费类型
						if i["free_type"]:
							# 区分有科室和没有科室
							if i["hosp_depart"]:
								requir_dict[i["business_type"]][i["report_type"]].setdefault((i["company"],i["prod_name"],i["hosp_depart"],i["order_type"],i["free_type"]), "")
								requir_dict[i["business_type"]][i["report_type"]][(i["company"],i["prod_name"],i["hosp_depart"],i["order_type"],i["free_type"])] = i["report_name"]
								brca_cnv_dict[i["business_type"]][i["report_type"]].setdefault((i["company"],i["prod_name"],i["hosp_depart"],i["order_type"],i["free_type"]), "")
								brca_cnv_dict[i["business_type"]][i["report_type"]][(i["company"],i["prod_name"],i["hosp_depart"],i["prod_name"],i["order_type"],i["free_type"])] = i["judge_brca_cnv"]
							else:
								requir_dict[i["business_type"]][i["report_type"]].setdefault((i["company"],i["prod_name"],i["order_type"],i["free_type"]), "")
								requir_dict[i["business_type"]][i["report_type"]][(i["company"],i["prod_name"],i["order_type"],i["free_type"])] = i["report_name"]
								brca_cnv_dict[i["business_type"]][i["report_type"]].setdefault((i["company"],i["prod_name"],i["order_type"],i["free_type"]), "")
								brca_cnv_dict[i["business_type"]][i["report_type"]][(i["company"],i["prod_name"],i["order_type"],i["free_type"])] = i["judge_brca_cnv"]
						# 不存在免费类型
						else:
							# 区分有科室和没有科室
							if i["hosp_depart"]:
								requir_dict[i["business_type"]][i["report_type"]].setdefault((i["company"],i["prod_name"],i["hosp_depart"],i["order_type"]), "")
								requir_dict[i["business_type"]][i["report_type"]][(i["company"],i["prod_name"],i["hosp_depart"],i["order_type"])] = i["report_name"]
								brca_cnv_dict[i["business_type"]][i["report_type"]].setdefault((i["company"],i["prod_name"],i["hosp_depart"],i["order_type"]), "")
								brca_cnv_dict[i["business_type"]][i["report_type"]][(i["company"],i["prod_name"],i["hosp_depart"],i["order_type"])] = i["judge_brca_cnv"]
							else:
								requir_dict[i["business_type"]][i["report_type"]].setdefault((i["company"],i["prod_name"],i["order_type"]), "")
								requir_dict[i["business_type"]][i["report_type"]][(i["company"],i["prod_name"],i["order_type"])] = i["report_name"]
								brca_cnv_dict[i["business_type"]][i["report_type"]].setdefault((i["company"],i["prod_name"],i["order_type"]), "")
								brca_cnv_dict[i["business_type"]][i["report_type"]][(i["company"],i["prod_name"],i["order_type"])] = i["judge_brca_cnv"]
					# 不存在order_type（进院）
					else:
						# 区分有科室和没有科室
						if i["hosp_depart"]:
							requir_dict[i["business_type"]][i["report_type"]].setdefault((i["company"],i["hosp_depart"],i["prod_name"]), "")
							requir_dict[i["business_type"]][i["report_type"]][(i["company"],i["hosp_depart"],i["prod_name"])] = i["report_name"]
							brca_cnv_dict[i["business_type"]][i["report_type"]].setdefault((i["company"],i["hosp_depart"],i["prod_name"]), "")
							brca_cnv_dict[i["business_type"]][i["report_type"]][(i["company"],i["hosp_depart"],i["prod_name"])] = i["judge_brca_cnv"]
						else:
							requir_dict[i["business_type"]][i["report_type"]].setdefault((i["company"],i["prod_name"]), "")
							requir_dict[i["business_type"]][i["report_type"]][(i["company"],i["prod_name"])] = i["report_name"]
							brca_cnv_dict[i["business_type"]][i["report_type"]].setdefault((i["company"],i["prod_name"]), "")
							brca_cnv_dict[i["business_type"]][i["report_type"]][(i["company"],i["prod_name"])] = i["judge_brca_cnv"]
			else:
					requir_dict[i["business_type"]][i["report_type"]].setdefault(i["prod_name"], "")
					requir_dict[i["business_type"]][i["report_type"]][i["prod_name"]] = i["report_name"]
					brca_cnv_dict[i["business_type"]][i["report_type"]].setdefault(i["prod_name"], "")
					brca_cnv_dict[i["business_type"]][i["report_type"]][i["prod_name"]] = i["judge_brca_cnv"]

	# 是否需要拼接基础模板对照表
	report_merge_dict = {}
	for i in Data:
		if str(int(i["status"])) == "0":
			if i["report_name"] not in report_merge_dict.keys():
				report_merge_dict.setdefault(i["report_name"],{})
			report_merge_dict[i["report_name"]] = {
				"merge_template" : str(int(i["merge_template"])),
				"temp_dynamic" : re.split(";", i["temp_dynamic"].replace("\n","").replace(" ","")) if i["temp_dynamic"] else "",
				"temp_fixed" : re.split(";", i["temp_fixed"].replace("\n","").replace(" ","")) if i["temp_fixed"] else ""
			}
	#print (requir_dict)

	return requir_dict, report_merge_dict, brca_cnv_dict
