#-*- coding:gbk -*-
import re
from libs.getRequirExcel import getRequirenment
from libs.getProdName_alias import alias_name
# 2026.01.21-新增适配临检LIMS来源订单
from libs.getRequirExcel import getRequirenment_v2
from libs.getConfig import get_department_name
# 2026.01.21-新增完成

'''
Discription
	
	匹配报告模板配置表，获取模板名称，用于报告填充。 

'''

def MatchReport(jsonDict, config):
	requir_dict, report_merge_dict, brca_cnv_dict = getRequirenment(config)
	i = jsonDict.get("sample_info")
	# 获取产品别名，转化为主要产品名，便于匹配模板
	alias_name_dict = alias_name(config)
	i["prod_names"] = alias_name_dict[i["prod_names"]] if i["prod_names"] in alias_name_dict.keys() else i["prod_names"]
	#匹配报告模板
	# 1. 临检
	 # 新增药企配置-匹配放前面-2023.11.17
	 # 1.0 先匹配"临检:定制:医院-产品-项目"
	 # 1.1 先匹配"临检:定制:医院-产品"
	 # 1.2 无定制的话汇总（2024.08.16新增-JY）匹配"临检:通用-简版-产品"， 否则匹配"临检:通用-完整-产品"
	# 2. 进院
	 # 2.1 先匹配定制"进院:定制:医院-科室-产品"，无的话匹配"进院:定制:医院-产品"
	 # 2.2 无定制的话匹配"进院:通用:产品"
	report_name = ""
	judge_brca_cnv = ""
	if i["report_module_type"] == "rummage":
		#print (i)
		# 20220627-这边做个兼容：临检样本需要匹配带后缀送检单位名称，而手动上传的订单缺失该项内容，若需要用手动上传订单出临检报告，会报错。
		# 这边兼容成缺少带后缀送检单位名称时，直接出临检完整版报告
		i["origin_company"] = i["origin_company"] if i["origin_company"] else i["company"]
		# 兼容完成-20220627
		# 新增药企配置-2023.11.17，单位名称先用company吧，现在还是手动上传的订单，不清楚对接LIMS后会是什么样
		if "product_name" in i.keys() and i["product_name"] and (i["company"], i["prod_names"], i["product_name"]) in requir_dict["rummage"]["CustomEdition"].keys():
			report_name = requir_dict["rummage"]["CustomEdition"].get((i["company"], i["prod_names"], i["product_name"]))
			judge_brca_cnv = brca_cnv_dict["rummage"]["CustomEdition"].get((i["company"], i["prod_names"], i["product_name"]))
		else:
		# 新增药企配置代码更新结束-2023.11.17
			#print ("yes")
			#print (requir_dict["rummage"]["CustomEdition"].keys())
			if (i["origin_company"], i["prod_names"]) in requir_dict["rummage"]["CustomEdition"].keys():
				#print ("yes")
				report_name = requir_dict["rummage"]["CustomEdition"].get((i["origin_company"], i["prod_names"]))
				judge_brca_cnv = brca_cnv_dict["rummage"]["CustomEdition"].get((i["origin_company"], i["prod_names"]))
			else:
				report_name = requir_dict["rummage"]["Universal_simple"].get(i["prod_names"], "") if re.search("汇总|JY", i["origin_company"]) else requir_dict["rummage"]["Universal_complete"].get(i["prod_names"], "")
				judge_brca_cnv = brca_cnv_dict["rummage"]["Universal_simple"].get(i["prod_names"], "") if re.search("汇总|JY", i["origin_company"]) else brca_cnv_dict["rummage"]["Universal_complete"].get(i["prod_names"], "")
	elif i["report_module_type"] == "hospital":
		if (i["company"], i["hosp_depart"], i["prod_names"]) in requir_dict["hospital"]["CustomEdition"].keys():
			report_name = requir_dict["hospital"]["CustomEdition"].get((i["company"], i["hosp_depart"], i["prod_names"]))
			judge_brca_cnv = brca_cnv_dict["hospital"]["CustomEdition"].get((i["company"], i["hosp_depart"], i["prod_names"]))
		elif (i["company"], i["prod_names"]) in requir_dict["hospital"]["CustomEdition"].keys():
			report_name = requir_dict["hospital"]["CustomEdition"].get((i["company"], i["prod_names"]))
			judge_brca_cnv = brca_cnv_dict["hospital"]["CustomEdition"].get((i["company"], i["prod_names"]))
		else:
			report_name = requir_dict["hospital"]["Universal"].get(i["prod_names"], "")
			judge_brca_cnv = brca_cnv_dict["hospital"]["Universal"].get(i["prod_names"], "")
	
	merge_template = report_merge_dict.get(report_name)

	return report_name, merge_template, judge_brca_cnv

def MatchReport_v2(jsonDict, config):
	'''
	Discription

	临检模板匹配规则更新
	所需字段：
		report_module_type：报告模块类型（临检：rummage，进院：hospital）
		company：送检单位名称
		prod_names：产品名称
		product_name：项目名称（可选）
		hosp_depart：医院科室
		order_type：订单类型，新增
		free_type：免费类型，新增
	规则：
		1. 临检：
			1.1 药企：
				匹配“临检：定制：医院-产品-项目”，有的话出药企模板	
			1.2 临检样本
				** 如果是临检项目-凑上机，则走进院模板匹配规则
				1.2.1 匹配“临检：定制：医院-产品-科室-free_type-order_type”， 有的话出医院科室免费类型定制模板
				1.2.2 匹配“临检：定制：医院-产品-free_type-order_type”，有的话出医院免费类型定制模板
				1.2.3 匹配“临检：定制：医院-产品-科室-order_type”，有的话出医院科室定制模板
				1.2.4 匹配“临检：定制：医院-产品-order_type”，有的话出医院定制模板
				1.2.5 上述未匹配到模板的，出通用规则模板
					临检项目：通用pdf
					汇总送检：通用pdf
					员工福利：通用pdf
					医学项目：通用pdf
					科研项目：通用pdf
					内部质控：通用pdf
					汇总送样：通用word（这个是汇总进院？）
					外部质控：通用word
		2. 进院：
			2.1 匹配“进院：定制：医院-产品-科室”，有的话出医院科室定制模板
			2.2 匹配“进院：定制：医院-产品”，有的话出医院定制模板
			2.3 上述未匹配到模板的，出进院通用模板
	'''
	requir_dict, report_merge_dict, brca_cnv_dict = getRequirenment_v2(config)
	i = jsonDict.get("sample_info")
	# 获取产品别名，转化为主要产品名，便于匹配模板
	alias_name_dict = alias_name(config)
	i["prod_names"] = alias_name_dict[i["prod_names"]] if i["prod_names"] in alias_name_dict.keys() else i["prod_names"]

	report_name = ""
	judge_brca_cnv = ""
	if i["report_module_type"] == "rummage":
		department_dict = get_department_name(config)
		i["institute_dept_id"] = i["institute_dept_id"] if "institute_dept_id" in i.keys() and i["institute_dept_id"] else ""
		i["hosp_depart"] = department_dict.get(i["institute_dept_id"], i["hosp_depart"]) 
		print ("标准科室ID：", i["institute_dept_id"], "标准科室名称：", department_dict.get(i["institute_dept_id"],""), "最终科室名：", i["hosp_depart"])
		# 药企-还是按之前的匹配规则
		if "product_name" in i.keys() and i["product_name"] and (i["company"], i["prod_names"], i["product_name"]) in requir_dict["rummage"]["CustomEdition"].keys():
			report_name = requir_dict["rummage"]["CustomEdition"].get((i["company"], i["prod_names"], i["product_name"]))
			judge_brca_cnv = brca_cnv_dict["rummage"]["CustomEdition"].get((i["company"], i["prod_names"], i["product_name"]))
		else:
			if i["free_type"]:
				# 如果是凑上机样本，当成进院的处理
				if i["order_type"] == "临检项目" and i["free_type"] == "凑上机":
					if (i["company"], i["hosp_depart"], i["prod_names"]) in requir_dict["hospital"]["CustomEdition"].keys():
						report_name = requir_dict["hospital"]["CustomEdition"].get((i["company"], i["hosp_depart"], i["prod_names"]))
						judge_brca_cnv = brca_cnv_dict["hospital"]["CustomEdition"].get((i["company"], i["hosp_depart"], i["prod_names"]))
					elif (i["company"], i["prod_names"]) in requir_dict["hospital"]["CustomEdition"].keys():
						report_name = requir_dict["hospital"]["CustomEdition"].get((i["company"], i["prod_names"]))
						judge_brca_cnv = brca_cnv_dict["hospital"]["CustomEdition"].get((i["company"], i["prod_names"]))
					else:
						report_name = requir_dict["hospital"]["Universal"].get(i["prod_names"], "")
						judge_brca_cnv = brca_cnv_dict["hospital"]["Universal"].get(i["prod_names"], "")
				# 其他
				else:
					if (i["company"], i["prod_names"], i["hosp_depart"], i["order_type"], i["free_type"]) in requir_dict["rummage"]["CustomEdition"].keys():
						report_name = requir_dict["rummage"]["CustomEdition"].get((i["company"], i["prod_names"], i["hosp_depart"], i["order_type"], i["free_type"]))
						judge_brca_cnv = brca_cnv_dict["rummage"]["CustomEdition"].get((i["company"], i["prod_names"], i["hosp_depart"], i["order_type"], i["free_type"]))
					elif (i["company"], i["prod_names"], i["order_type"], i["free_type"]) in requir_dict["rummage"]["CustomEdition"].keys():
						report_name = requir_dict["rummage"]["CustomEdition"].get((i["company"], i["prod_names"], i["hosp_depart"], i["order_type"], i["free_type"]))
						judge_brca_cnv = brca_cnv_dict["rummage"]["CustomEdition"].get((i["company"], i["prod_names"], i["hosp_depart"], i["order_type"], i["free_type"]))
			else:
				if (i["company"], i["prod_names"], i["hosp_depart"], i["order_type"]) in requir_dict["rummage"]["CustomEdition"].keys():
					report_name = requir_dict["rummage"]["CustomEdition"].get((i["company"], i["prod_names"], i["hosp_depart"], i["order_type"]))
					judge_brca_cnv = brca_cnv_dict["rummage"]["CustomEdition"].get((i["company"], i["prod_names"], i["hosp_depart"], i["order_type"]))
				elif (i["company"], i["prod_names"], i["order_type"]) in requir_dict["rummage"]["CustomEdition"].keys():
					report_name = requir_dict["rummage"]["CustomEdition"].get((i["company"], i["prod_names"], i["order_type"]))
					judge_brca_cnv = brca_cnv_dict["rummage"]["CustomEdition"].get((i["company"], i["prod_names"], i["order_type"]))
			if not report_name:
				# 默认模板
				# 2026.01.26-汇总送样改为汇总送检
				# 2026.03.19-汇总送检默认出pdf
				if i["order_type"] in ["外部质控", "汇总进院"]:
					report_name = requir_dict["rummage"]["Universal_simple"].get(i["prod_names"], "")
					judge_brca_cnv = brca_cnv_dict["rummage"]["Universal_simple"].get(i["prod_names"], "")
				else:
					report_name = requir_dict["rummage"]["Universal_complete"].get(i["prod_names"], "")
					judge_brca_cnv = brca_cnv_dict["rummage"]["Universal_complete"].get(i["prod_names"], "")
	elif i["report_module_type"] == "hospital":
		if (i["company"], i["hosp_depart"], i["prod_names"]) in requir_dict["hospital"]["CustomEdition"].keys():
			report_name = requir_dict["hospital"]["CustomEdition"].get((i["company"], i["hosp_depart"], i["prod_names"]))
			judge_brca_cnv = brca_cnv_dict["hospital"]["CustomEdition"].get((i["company"], i["hosp_depart"], i["prod_names"]))
		elif (i["company"], i["prod_names"]) in requir_dict["hospital"]["CustomEdition"].keys():
			report_name = requir_dict["hospital"]["CustomEdition"].get((i["company"], i["prod_names"]))
			judge_brca_cnv = brca_cnv_dict["hospital"]["CustomEdition"].get((i["company"], i["prod_names"]))
		else:
			report_name = requir_dict["hospital"]["Universal"].get(i["prod_names"], "")
			judge_brca_cnv = brca_cnv_dict["hospital"]["Universal"].get(i["prod_names"], "")
	
	merge_template = report_merge_dict.get(report_name)

	return report_name, merge_template, judge_brca_cnv

# 判断使用哪套模板匹配规则
def choose_template(jsonDict, config):
	'''
	根据jsonDict中的信息，判断使用哪套模板匹配规则
	'''
	i = jsonDict.get("sample_info")
	if i["report_module_type"] == "rummage" and "order_type" in i.keys() and i["order_type"]:
		print ("#---------存在业务类型字段，使用LIMS模板匹配规则---------#")
		report_name, merge_template, judge_brca_cnv = MatchReport_v2(jsonDict, config)
	else:
		print ("#---------不存在业务类型字段，使用BIMS模板匹配规则---------#")
		report_name, merge_template, judge_brca_cnv = MatchReport(jsonDict, config)
	return report_name, merge_template, judge_brca_cnv