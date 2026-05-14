#-*- coding:gbk -*-
import re
import os
import json
import datetime
import argparse
import time
import copy
import lxml
from docxtpl import DocxTemplate, InlineImage, R
from docx.shared import Mm
from bin import getSampleInfo
#from bin import getVar_BRCA
from bin import getQC
from bin import getTemplate
from bin import getApprovalDrug
from bin import getClinicTrial
from bin import getPDL1
from bin import getMSI
from bin import getTMB
from bin import getVar
from bin import getRNAexp
from bin import getChemo
from bin import getRefence
from bin import getGEP
from bin import getHRD
from bin import getTME
#from bin import getVar_HRR_SHSY
from bin import getApprovalRegimen
from bin import baseTemplate
from libs.getProdName_alias import alias_name
from bin import getImage
import customize_filters
import io
from bin import removeMSSEvi
from bin import removeTP53Evi
from bin import removeEvi
# 新增配置-2025.05.06
# 2025.07.23-新增get_international_department_risk
from libs.getConfig import get_risk_inter, get_risk_table, get_risk_management, get_hx_breast_templage_significance, get_international_department_risk
# 2025.06.03新增
from bin import removeKRASEvi
from merge_product import zszl_brcav1_tbrca, zszl_brcav1_hrd
# 2025.10.23新增
from bin import get_summary_for_specialcomany
# 2025.12.02-新增MRD结果
from bin import getMRD
# 2026.01.09-新增华大HD处理
from merge_product import zszl_brcav1_huada
# 2026.01.21-临检LIMS来源样本捕获项目文库总量使用dna_pre_library_qty/rna_pre_library_qty，临时方案直接赋值给library_qty
from  bin  import libQC_stran
# 2026.01.21-新增完成
# 2026.02.04-西安交大一-新增brca_v1+bhd合并报告
from merge_product import xajdy_brcav1_bhd
# 2026.02.26-山东齐鲁-新增HRD+gBRCA合并报告
from merge_product import sdql_hrd_gbrca
# 2026.04.03-吉大一MP DNA + RNAseq
from merge_product import jdyy_masterdna_rnaseq

def get_data(json_name, outfile, config, report_template, outjson, image):
	print ("开始生成填充数据！", datetime.datetime.now())
	data = {}
	# 解析json文件
	with open(outfile+"/"+json_name+".json", "r", encoding='utf-8') as file_json:
		jsonDict = json.load(file_json)
		jsonDict["sample_info"]["raw_prod_names"] = jsonDict["sample_info"]["prod_names"]
	
	# 2025.09.29-如果merge_order不为空，判定为多项目出一个报告，格式需要特殊处理
	# 2025.12.02-更新规则，MRD需要填写身份证/接收时间，可能导致模板匹配判定出错，这边兼容下
	report_name_judge_merge = ""
	if "merge_order" in jsonDict.keys() and jsonDict["merge_order"] and len(jsonDict["merge_order"]["cnv"]["sample_id_list"]) >= 2:
		merge_sample_info = {}
		merge_product = []
		for sample in jsonDict["merge_order"]["cnv"]["sample_id_list"]:
			# 2025.12.02-判断是否存在json文件
			if os.path.exists(outfile+"/"+sample["order_code"]+".json"):
				with open(outfile+"/"+sample["order_code"]+".json", "r", encoding='utf-8') as sample_json:
					sample_jsonDict = json.load(sample_json)
					merge_sample_info = sample_jsonDict["sample_info"]
					merge_product.append(sample_jsonDict["sample_info"]["prod_names"])
		merge_sample_info["prod_names"] = "，".join(sorted(merge_product))
		# 2026.01.21-兼容BIMS/LIMS来源订单模板匹配
		#report_name, merge_template, judge_brca_cnv = getTemplate.MatchReport({"sample_info" : merge_sample_info}, config)
		report_name, merge_template, judge_brca_cnv = getTemplate.choose_template({"sample_info" : merge_sample_info}, config)
		# 2026.01.21-兼容完成
		if set(merge_product) == set(["BRCA1/2（扩增子）", "BRCA1/BRCA2（组织）"]) and merge_sample_info["company"] == "中山大学附属肿瘤医院":
			data = zszl_brcav1_tbrca.MergeReport_main(jsonDict["merge_order"]["cnv"]["sample_id_list"], outfile)
		if set(merge_product) == set(["BRCA1/2（扩增子）", "HRD Complete（组织）"]) and merge_sample_info["company"] == "中山大学附属肿瘤医院":
			data = zszl_brcav1_hrd.MergeReport_main(jsonDict["merge_order"]["cnv"]["sample_id_list"], outfile)
		if set(merge_product) == set(["BRCA1/BRCA2（全血）", "BRCA1/BRCA2（组织 全血）"]) and merge_sample_info["company"] in ["西安交通大学第一附属医院", "西安交通大学医学院第一附属医院"]:
			data = xajdy_brcav1_bhd.MergeReport_main(jsonDict["merge_order"]["cnv"]["sample_id_list"], outfile)
		if set(merge_product) == set(["BRCA1/BRCA2（全血）", "HRD Complete（组织）"]) and merge_sample_info["company"] in ["山东大学齐鲁医院"]:
			data = sdql_hrd_gbrca.MergeReport_main(jsonDict["merge_order"]["cnv"]["sample_id_list"], outfile)

		report_name_judge_merge = report_name


	# 2026.04.03-吉大一MP DNA + RNAseq，merge_order.cnv.sample_id_list无法返回RNAseq的id，所以只能通过sample.relation_order_code来处理，仅考虑BIMS来源，LIMS等订单类型确认了再加
	if jsonDict["sample_info"]["report_module_type"] == "clinical" and jsonDict["sample_info"]["origin_company"] == "吉林大学第一医院-SF-肺癌":
		if "relation_order_code" in jsonDict["sample_info"].keys() and jsonDict["sample_info"]["relation_order_code"]:
			merge_sample_info = jsonDict["sample_info"]
			merge_product = []
			merge_product.append(jsonDict["sample_info"]["prod_names"])
			#jdyy_merge_prod_list = []
			#jdyy_merge_prod_list.append(jsonDict["sample_info"]["prod_names"])
			relation_jsonDict = {}
			if os.path.exists(outfile+"/"+jsonDict["sample_info"]["relation_order_code"][0]+".json"):
				with open(outfile+"/"+jsonDict["sample_info"]["relation_order_code"][0]+".json", "r", encoding='utf-8') as sample_json:
					relation_jsonDict = json.load(sample_json)
					relation_jsonDict["sample_info"]["raw_prod_names"] = relation_jsonDict["sample_info"]["prod_names"]
					merge_product.append(relation_jsonDict["sample_info"]["prod_names"])
			merge_sample_info["prod_names"] = "，".join(sorted(merge_product))
			merge_sample_info["report_module_type"] = "rummage" if merge_sample_info["report_module_type"] == "clinical" else merge_sample_info["report_module_type"]
			report_name, merge_template, judge_brca_cnv = getTemplate.choose_template({"sample_info" : merge_sample_info}, config)
			if report_name and set(merge_product) == set(["Master Panel（组织）", "RNASeq（肉瘤）"]):
				data = jdyy_masterdna_rnaseq.MergeReport_main(relation_jsonDict, jsonDict, report_name)
				report_name_judge_merge = report_name
	# 2026.04.03-新增完成


	# 2025.12.02-改为判断是否匹配到report_name，未匹配到的话按照单项目出报告
	if not report_name_judge_merge:
	# 2025.09.29-兼容完成
	
		# 系统将rummage改为clinical了，为了报告脚本代码改动最小化，在这边做转化
		jsonDict["sample_info"]["report_module_type"] = "rummage" if jsonDict["sample_info"]["report_module_type"] == "clinical" else jsonDict["sample_info"]["report_module_type"]
		# 将产品别名进行转换
		alias_name_dict = alias_name(config)
		jsonDict["sample_info"]["prod_names"] = alias_name_dict[jsonDict["sample_info"]["prod_names"]] if jsonDict["sample_info"]["prod_names"] in alias_name_dict.keys() else jsonDict["sample_info"]["prod_names"]
		# v4删除了tumor_names_cn，为了报告生成不报错，这边加兼容-2024.03.26
		jsonDict["sample_info"]["tumor_names_cn"] = jsonDict["sample_info"]["tumor_type"]
		# 这边做个MSS开关，T时需要删除MSS的证据及治疗方案介绍中的MSS-2024.05.07
		prod_list = ["BPTM Plus（组织）"]
		remove_MSS = "T"
		if remove_MSS == "T" and jsonDict["sample_info"]["prod_names"] in prod_list:
			jsonDict = removeMSSEvi.remove_MSS_evi(jsonDict)
		# 如果符合西安交大一相关条件时，过滤TP53 EGFR-TKIs证据-2025.03.19
		# 加个开关，T时删除TP53相关证据
		remove_TP53 = "T"
		if remove_TP53 == "T":
			jsonDict = removeTP53Evi.remove_TP53_evi(jsonDict)
		# 增加完成-2025.03.19
		# 符合西安交大一相关条件时，KRAS EGFR-TKIs证据降级，A级预后证据不参与评级-2025.06.03
		remove_KRAS = "T"
		if remove_KRAS:
			jsonDict = removeKRASEvi.remove_KRAS_evi(jsonDict)
		# 2025.06.12-新增，重庆西南116 TP53删除迈华替尼
		remove_cqxn_TP53 = "T"
		if remove_cqxn_TP53:
			jsonDict = removeTP53Evi.cqzn_remove_TP53_evi(jsonDict)
		# 2025.06.12-新增完成
		# 2025.08.20-新增，重庆西南116 PIK3CA非获批位点不展示证据
		remove_cqxn_PIK3CA = "T"
		if remove_cqxn_PIK3CA:
			jsonDict = removeEvi.cqxn_remove_PIK3CA_evi(jsonDict)
		# 2025.08.20-新增完成
		# 2026.02.05-新增，西安交大一体系模板，删除KRAS、NRAS中的呋喹替尼、瑞戈替尼、贝伐珠单抗（引用机构需要时CSCO）
		remove_xajdy_ras_evi = "T"
		if remove_xajdy_ras_evi:
			jsonDict = removeEvi.xajdy_remove_ras_evi(jsonDict)
		# 2026.02.05-新增完成
		# 2026.04.01-重庆西南删除CDKN2A Inactivating Mutation中的奥希替尼证据
		remove_cqxn_CDKN2A_evi = "T"
		if remove_cqxn_CDKN2A_evi:
			jsonDict = removeEvi.cqxn_remove_CDKN2A_evi(jsonDict)
		# 2026.04.01-新增完成

		# 模板选择
		# 2026.01.21-兼容BIMS/LIMS来源订单模板匹配
		#report_name, merge_template, judge_brca_cnv = getTemplate.MatchReport(jsonDict, config)
		report_name, merge_template, judge_brca_cnv = getTemplate.choose_template(jsonDict, config)
		# 2026.01.21-兼容完成

		# 没有设置的默认用MLPA
		print (judge_brca_cnv)
		data["judge_brca_cnv"] = judge_brca_cnv if judge_brca_cnv else "mlpa"
		#print (report_name)
		#print (merge_template)
		data["sample"] = getSampleInfo.getSample(jsonDict)
		data["sample"]["report_name"] = report_name
		data["qc"], data["lib_quality_control"] = getQC.getJsonQC(jsonDict)
		# 2026.01.21-临检LIMS来源捕获项目文库总量使用dna_pre_library_qty/rna_pre_library_qty，临时方案直接赋值给library_qty
		if "order_type" in data["sample"].keys() and data["sample"]["order_type"] and data["sample"]["report_module_type"] == "rummage":
			data["lib_quality_control"] = libQC_stran.clinical_libraty_qty(jsonDict["sample_info"]["report_module_type"], jsonDict["sample_info"]["prod_names"], data["lib_quality_control"])
		# 2026.01.21-新增完成
		# ---------------------------------------------------------------------------------------------------
		# 2026.01.29-LIMS肿瘤细胞含量返回规则有改动
		# 玻片病理结果无法传到订单中，改为从湿实验质控里获取，蜡卷还是从订单里获取
		# 为了不影响报告生成，优先从湿实验质控里获取（lib传递给sample，模板就不用改动了），再从订单里获取
		if "order_type" in data["sample"].keys() and data["sample"]["order_type"] and data["sample"]["report_module_type"] == "rummage":
			if "lib_dna_qc" in data["lib_quality_control"].keys() and \
				data["lib_quality_control"]["lib_dna_qc"] and \
				"tumor_content" in data["lib_quality_control"]["lib_dna_qc"].keys() and \
				data["lib_quality_control"]["lib_dna_qc"]["tumor_content"]:

				data["sample"]["tumor_content"] = data["lib_quality_control"]["lib_dna_qc"]["tumor_content"]
				data["sample"]["tumor_content_num"] = data["sample"]["tumor_content"].replace("%", "") if \
													 "tumor_content" in data["sample"].keys() and \
													  data["sample"]["tumor_content"] else \
													  ""
				print ("LIMS临检订单，肿瘤细胞含量，湿实验：{0}，最终：{1}，num值：{2}".format(
				data["lib_quality_control"]["lib_dna_qc"]["tumor_content"],
				data["sample"]["tumor_content"],
				data["sample"]["tumor_content_num"]
				))

			if "lib_dna_qc" in data["lib_quality_control"].keys() and \
				data["lib_quality_control"]["lib_dna_qc"] and \
				"tumor_content_macrodissection_performed" in data["lib_quality_control"]["lib_dna_qc"].keys() and \
				data["lib_quality_control"]["lib_dna_qc"]["tumor_content_macrodissection_performed"]:

				data["sample"]["tumor_content_macrodissection_performed"] = data["lib_quality_control"]["lib_dna_qc"]["tumor_content_macrodissection_performed"]
				data["sample"]["tumor_content_macrodissection_performed_num"] = data["sample"]["tumor_content_macrodissection_performed"].replace("%", "") if \
																				"tumor_content_macrodissection_performed" in data["sample"].keys() and \
																				data["sample"]["tumor_content_macrodissection_performed"] else \
																				""
				print ("LIMS临检订单，富集后肿瘤细胞含量，湿实验：{0}，最终：{1}，num值：{2}".format(
				data["lib_quality_control"]["lib_dna_qc"]["tumor_content_macrodissection_performed"],
				data["sample"]["tumor_content_macrodissection_performed"],
				data["sample"]["tumor_content_macrodissection_performed_num"]
				))
			# 体系项目的有富集后看富集后肿瘤细胞含量，没有看原始肿瘤细胞含量，胚系项目保留两个肿瘤细胞含量
			somatic_prod_list = ["10基因（组织）", "CRC12-MSI", "Pan116（组织）", "Master Panel（组织）", "Classic Panel", "TC21（组织）", \
					  	 "GA18（组织）", "LC76（组织）", "CRC25（组织）", "Master Panel（组织-综合）", "Classic Panel-胃", "Classic Panel-胃肠", \
						 "Classic Panel-甲状腺", "Classic Panel-肝胆", "Classic Panel-骨", "Classic Panel-黑色素", "Classic Panel-头颈", \
						 "Pan116（组织）-胰腺", "Pan116（组织）-乳腺", "OncoPro（组织）", "Master Panel V2（组织）", "Master Panel V1（组织）", \
						 "RNASeq（肉瘤）", "Pan116（组织）-黑色素瘤", "Pan116（组织）-头颈", "OncoPro（组织）-妇科肿瘤", "Master Panel V1（组织_北医）", \
						 "Classic Panel 200（组织）"]
			if data["sample"]["prod_names"] in somatic_prod_list:
				data["sample"]["tumor_content"] = data["sample"]["tumor_content_macrodissection_performed"] if \
												  "tumor_content_macrodissection_performed" in data["sample"].keys() and \
												  data["sample"]["tumor_content_macrodissection_performed"] else \
												  data["sample"]["tumor_content"]
				data["sample"]["tumor_content_num"] = data["sample"]["tumor_content"].replace("%", "") if \
													 "tumor_content" in data["sample"].keys() and \
													  data["sample"]["tumor_content"] else \
													  ""
		# 2026.01.29-新增完成
		# ---------------------------------------------------------------------------------------------------------------
		data["drug"] = getApprovalDrug.getDrug(jsonDict)
		data["therapeutic_regimen"] = getApprovalRegimen.getRegimen(jsonDict)
		data["clinic_trial"] = getClinicTrial.getClinic(jsonDict)
		data["pdl1"] = getPDL1.getPDL1(jsonDict)
		data["msi"] = getMSI.getMSI(jsonDict, config)
		data["tmb"] = getTMB.getTMB(jsonDict, config)
		data["var"], data["var_brca"], data["var_hrr_shsy"] = getVar.getVar(jsonDict, config, report_name)
		data["rna_exp"] = getRNAexp.getRNA_exp(jsonDict)
		data["chemo"] = getChemo.getchemo(jsonDict, config)
		data["gep"] = getGEP.getGEP(jsonDict)
		data["hrd"] = getHRD.getHRD(jsonDict, data["var"]["ec_type"]["BRCA1_level12"]+data["var"]["ec_type"]["BRCA2_level12"], config)
		data["tme"], data["tme_score"] = getTME.getTME(jsonDict)
		# 用于浙肿BRCA判定是否有检出CNV信号
		data["judge_CNV"] = "T" if jsonDict["cnv"] else ""
		# 参考文献
		data["refer"] = {}
		data["refer"]["fixed"] = getRefence.getfixed_refer(report_name, data["sample"]["tumor_list"], data["var"]["knb"], data["msi"], config)
		data["refer"]["dynamic"] = getRefence.getdynamic_refer(jsonDict, data["var"], data["msi"], data["hrd"], data["var_brca"])
		# 新增生信输出的所有化疗结果-2024.07.31
		data["PGx_analysis_result"] = jsonDict["PGx_analysis_result"] if "PGx_analysis_result" in jsonDict.keys() and jsonDict["PGx_analysis_result"] else []
		# 新增hd结果-2024.10.12
		data["hd"] = jsonDict["hd"] if "hd" in jsonDict.keys() and jsonDict["hd"] else []
		# 2025.05.06-新增配置
		data["hx_150_gene_risk_inter"] = get_risk_inter(config)
		data["hx_150_gene_risk_table"] = get_risk_table(config)
		data["hx_150_gene_risk_management"] = get_risk_management(config)
		data["hx_breast_template_significance"] = get_hx_breast_templage_significance(config)
		# 新增结束
		# 2025.07.23-新增配置
		data["international_department_risk"] = get_international_department_risk(config)
		# 2025.07.23-新增完成
		# 增加返回refer原始结果-2025.05.16
		data["refer_original"] = jsonDict["refer"]
		# 新增完成-2025.04.29
		# 2025.09.24-新增report_audit
		data["report_audit"] = jsonDict["report_audit"] if "report_audit" in jsonDict.keys() else {}
		# 2025.09.24-新增完成
		# 2025.12.02-新增MRD结果
		data["mrd"] = getMRD.getMRD(jsonDict)
		# 2025.12.02-新增完成
		# 2026.01.09-新增华大
		if "hua_results" in jsonDict.keys() and jsonDict["hua_results"]:
			data["huada_data"] = zszl_brcav1_huada.get_huada_data(jsonDict["hua_results"])
		# 2026.01.09-新增完成


	# 加一个PQCC hrd_txt 临时使用-2024.02.04
	#data["hrd_txt"] = jsonDict["hrd_txt"] if "hrd_txt" in jsonDict.keys() else {}
	# 2024.02.04-用完删掉

	### 厦门市一输出refre到json-2023.07.17
	if re.search("XMSY", report_name):
		refer_forXMSY = []
		# 1. CP40
		if data["sample"]["prod_names"] == "Classic Panel":
			refer_forXMSY.append("美国国家综合癌症网络（NCCN\u00AE） 肿瘤临床实践指南。")
			for i in data["refer"]["fixed"]+data["refer"]["dynamic"]["s_var12"]+data["refer"]["dynamic"]["msi"]+data["refer"]["dynamic"]["knb"]:
				if i not in refer_forXMSY:
					refer_forXMSY.append(i)
		# 2. BRCA
		elif data["sample"]["prod_names"] in ["BRCA1/BRCA2（全血）", "BRCA1/BRCA2（组织）", "BRCA1/BRCA2（组织 全血）"]:
			#brca_fix_refer = ["National Comprehensive Cancer Network (2020) NCCN Clinical Practice Guidelines in Oncology-Genetic/Familial High-Risk Assessment: Breast, Ovarian, and Pancreatic. (version 1.2020)[EB/OL]. http://www.nccn.org",
			#				  "National Comprehensive Cancer Network (2020) NCCN Clinical Practice Guidelines in Oncology-Breast Cancer Risk Reduction. (version 1.2020)[EB/OL]. http://www.nccn.org",
			#				  "National Comprehensive Cancer Network (2020) NCCN Clinical Practice Guidelines in Oncology-Ovarian Cancer Including Fallopian Tube Cancer and Primary Peritoneal Cancer. (version 1.2020)[EB/OL]. http://www.nccn.org",
			#				  "National Comprehensive Cancer Network (2020) NCCN Clinical Practice Guidelines in Oncology-Breast Cancer. (version 4.2020)[EB/OL]. http://www.nccn.org",
			#				  "National Comprehensive Cancer Network (2020) NCCN Clinical Practice Guidelines in Oncology-Prostate Cancer. (version 2.2020)[EB/OL]. http://www.nccn.org",
			#				  "National Comprehensive Cancer Network (2020) NCCN Clinical Practice Guidelines in Oncology-Pancreatic Adenocarcinoma. (version 1.2020)[EB/OL]. http://www.nccn.org",
			#				  "中国医师协会精准治疗委员会乳腺癌专业委员会, 中华医学会肿瘤学分会乳腺肿瘤学组, 中国抗癌协会乳腺癌专业委员会, 等. (2018) 中国乳腺癌患者BRCA1/2基因检测与临床应用专家共识(2018年版). 中国癌症杂志 28(10):787-800."]
			# BRCA参考文献更新-2024.05.17
			brca_fix_refer = [
				"中华医学会病理学分会, 国家病理质控中心. BRCA1/2数据解读中国专家共识（2021版）[J] .中华病理学杂志, 2021, 50(6): 565-571.",
				"中国医师协会精准治疗委员会乳腺癌专业委员会, 中华医学会肿瘤学分会乳腺肿瘤学组, 中国抗癌协会乳腺癌专业委员会. 中国乳腺癌患者BRCA1/2基因检测与临床应用专家共识（2018年版）[J].中国癌症杂志, 2018, 28(10): 787-800.",
				"中国抗癌协会肿瘤标志专业委员会遗传性肿瘤标志物协作组, 中华医学会病理学分会分子病理学组. 同源重组修复缺陷临床检测与应用专家共识（2021版）[J].中国癌症防治杂志, 2021, 13(4): 329-338.",
				"中华医学会妇科肿瘤学分会. 卵巢癌PARP抑制剂临床应用指南（2022版）[J].现代妇产科进展,2022,31(8):561-572.",
				"中国抗癌协会泌尿男生殖系肿瘤专业委员会, 中国临床肿瘤学会前列腺癌专业委员会. 中国前列腺癌患者基因检测专家共识（2020年版）[J].中国癌症杂志, 2020, 30(7): 551-560.",
				"NCCN Clinical Practice Guidelines in Oncology-Genetic/Familial High-Risk Assessment: Breast, Ovarian, and Pancreatic. (version 2.2024)[EB/OL]. http://www.nccn.org",
				"NCCN Clinical Practice Guidelines in Oncology-Breast Cancer Risk Reduction. (version 1.2024)[EB/OL]. http://www.nccn.org",
				"NCCN Clinical Practice Guidelines in Oncology-Ovarian Cancer Including Fallopian Tube Cancer and Primary Peritoneal Cancer. (version 2.2023)[EB/OL]. http://www.nccn.org",
				"NCCN Clinical Practice Guidelines in Oncology-Breast Cancer. (version 4.2023)[EB/OL]. http://www.nccn.org",
				"NCCN Clinical Practice Guidelines in Oncology-Prostate Cancer. (version 4.2023)[EB/OL]. http://www.nccn.org",
				"NCCN Clinical Practice Guidelines in Oncology-Pancreatic Adenocarcinoma. (version 2.2023)[EB/OL]. http://www.nccn.org",
				"Richards S, Aziz N, Bale S, et al. (2015) Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology. Genet Med.;17(5):405-424. [PMID: 25741868]",
				"王秋菊, 沈亦平, 邬玲仟, 陈少科, 陈子江, 方向东, 傅松滨, 龚瑶琴, 黄国英, 黄国宁, 黄荷凤, 黄山, 郝晓柯 12, 冀小平, 李红, 梁波, 廖灿, 乔杰, 苏海翔, 魏军, 王磊, 王树玉, 王晓红, 邢清和, 徐湘民, 袁慧军, 杨正林, 周从容, 周文浩, 曾勇, 张学军, 黄涛生, 郑茜, 秦胜营, 于世辉, 关静, 王洪阳, 王大勇, 赵立东, 王慧君, 孔令印, 宣黎明, 冒燕, 祝轶君, 徐君玲, 王剑青, 王莉, 赵婷, 秦一丁, 夏滢颖, 樊丽霞, 赵丁丁, 邱浩, 贺林.  遗传变异分类标准与指南[J]. 中国科学：生命科学, 2017, 47(6): 668-688.",
				"Li MM, Datto M, Duncavage EJ, Kulkarni S, et al. (2017) Standards and Guidelines for the Interpretation and Reporting of Sequence Variants in Cancer: A Joint Consensus Recommendation of the Association for Molecular Pathology American Society of Clinical Oncology, and College of American Pathologists. J Mol Diagn. Jan;19(1):4-23. [PMID: 27993330]",
				"二代测序临床报告解读专家组. 二代测序临床报告解读指引[J]. 循证医学, 2020, 20(4): 193-202."
			]
			refer_forXMSY.extend(brca_fix_refer)
			for i in data["var_brca"]["refer"]["var"]+data["var_brca"]["refer"]["gene"]+data["var_brca"]["refer"]["evi"]:
				if i not in refer_forXMSY:
					refer_forXMSY.append(i)
		# 3. LC10
		elif data["sample"]["prod_names"] in ["10基因（组织）", "10基因（血液）"]:
			lc10_fix_refer = ["美国国家综合癌症网络（NCCN\u00AE） 肿瘤临床实践指南。",
							  "Debyani Chakravarty, Jianjiong Gao, Sarah M Phillips et al.(2017) OncoKB: A Precision Oncology Knowledge Base. JCO Precis Oncol. 2017 Jul;2017:PO.17.00011.[PMID: 28890946]",
							  "Marilyn M Li, Michael Datto, Eric J Duncavage.(2017) Standards and Guidelines for the Interpretation and Reporting of Sequence Variants in Cancer: A Joint Consensus Recommendation of the Association for Molecular Pathology, American Society of Clinical Oncology, and College of American Pathologists. J Mol Diagn.2017 Jan;19(1):4-23.[PMID: 27993330]"]
			refer_forXMSY.extend(lc10_fix_refer)
			for i in data["refer"]["fixed"]+data["refer"]["dynamic"]["s_var12"]+data["refer"]["dynamic"]["knb"]:
				if i not in refer_forXMSY:
					refer_forXMSY.append(i)
		# 4. PAN116
		elif data["sample"]["prod_names"] in ["Pan116（组织）", "Pan116（血液）"]:
			refer_forXMSY.append("美国国家综合癌症网络（NCCN\u00AE） 肿瘤临床实践指南。")
			for i in data["refer"]["fixed"]+data["refer"]["dynamic"]["s_var12"]+data["refer"]["dynamic"]["g_var45"]+data["refer"]["dynamic"]["knb"]+data["refer"]["dynamic"]["msi"]:
				if i not in refer_forXMSY:
					refer_forXMSY.append(i)
		refer_dataJson = json.dumps({"refer" : refer_forXMSY}, ensure_ascii = False)
		with open(outfile+"/"+json_name+"_refer.json", "w", encoding = "utf-8") as outFile:
			outFile.write(refer_dataJson)
	### 厦门市一输出refre到json结束-2023.07.17

	### 厦门市一CP200 （无定制模板，目前用的是临检通用版）输出refer到json-2025.11.07
	# 产品名 OncoPro（组织）、Classic Panel 200（组织），进院
	if data["sample"]["prod_names"] in ["OncoPro（组织）", "Classic Panel 200（组织）"] and data["sample"]["company"] == "厦门大学附属第一医院" and data["sample"]["report_module_type"] == "hospital":
		refer_forXMSY = []
		refer_forXMSY.append("美国国家综合癌症网络（NCCN） 肿瘤临床实践指南")
		for i in data["refer"]["fixed"]+data["refer"]["dynamic"]["s_var12"]+data["refer"]["dynamic"]["s_var_onco_nodrug"]+data["refer"]["dynamic"]["msi"] +data["refer"]["dynamic"]["knb"]:
			if i not in refer_forXMSY:
				refer_forXMSY.append(i)
		refer_dataJson = json.dumps({"refer" : refer_forXMSY}, ensure_ascii = False)
		with open(outfile+"/"+json_name+"_refer.json", "w", encoding = "utf-8") as outFile:
			outFile.write(refer_dataJson)
	### 厦门市一CP200 （无定制模板，目前用的是临检通用版）输出refer到json结束-2025.11.07

	### 福建肿瘤LC10输出结果解读到json（仅进院）-2025.07.01
	if data["sample"]["prod_names"] == "10基因（组织）" and re.search("FJZL", report_name) and data["sample"]["report_module_type"] == "hospital":
		inter_forFJZL = []
		if data["var"]["summary_result"]["format9_forLC10_FJZL"]["var"]:
			for i in data["var"]["summary_result"]["format9_forLC10_FJZL"]["var"]:
				inter_forFJZL.append(i["gene_info"])
				inter_forFJZL.append(i["rpt_info"])
				inter_forFJZL.append(i["var_info"])
		if data["var"]["knb"]:
			inter_forFJZL.append("KRAS属于原癌基因，编码的蛋白在细胞内的信号通路中发挥信号转导作用。KRAS突变主要发生在Exon2或3的第12、13和61号密码子上，这些位点位于蛋白的GTP酶结构域。KRAS突变在多种肿瘤中均有发生，包括肺癌，结直肠癌和胰腺癌等等（PMID: 18794081, PMID: 19679400, PMID: 20952405）。NRAS基因属于原癌基因，编码的RAS蛋白在多种细胞信号通路中起着信号转导作用，在细胞的生存与增殖等活动中发挥重要作用。目前，已经在多种肿瘤中发现NRAS基因突变，如黑色素瘤、结直肠癌、甲状腺癌等。BRAF，丝氨酸/苏氨酸蛋白激酶B-raf，为Raf丝氨酸/苏氨酸蛋白激酶家族的成员之一，通过MAP激酶途径发出信号以调节细胞增殖和细胞生长（PMID: 24737949, PMID: 29540830）。已在多种癌症中鉴定出BRAF突变，包括结直肠癌（PMID: 30122982），肺癌（PMID: 29729495），甲状腺（PMID: 12970315）和黑色素瘤（PMID: 24737949），并且还发现了许多突变证明会产生耐药性（PMID: 27478040）。")
			inter_forFJZL.append("NCCN结直肠癌指南中认为，所有转移性结直肠癌患者均应进行RAS基因（KRAS和NRAS）及BRAF基因突变检测；并指出，KRAS、NRAS基因突变（2/3/4号外显子中）的结直肠癌患者不能用西妥昔单抗或帕尼单抗进行治疗；而BRAF V600E突变的患者除非同时使用BRAF抑制剂，否则极有可能对西妥昔单抗和帕尼单抗无反应。CSCO结直肠癌指南（2022）中，对RAS和BRAF均野生型的患者，推荐使用西妥昔单抗和贝伐珠单抗联合不同化疗；对RAS或BRAF突变的患者，推荐贝伐珠单抗联合不同化疗；对RAS野生、BRAF V600E突变的患者，推荐伊立替康+西妥昔单抗+维莫非尼、BRAF抑制剂+西妥昔单抗±MEK抑制剂。")
		fjzl_knb_regimen = []
		if data["var"]["knb"]:
			if data["var"]["knb"]["evi_sum"]["regimen_S_str"]:
				fjzl_knb_regimen.append("对{0}敏感".format(data["var"]["knb"]["evi_sum"]["regimen_S_str"]))
			if data["var"]["knb"]["evi_sum"]["regimen_R_str"]:
				fjzl_knb_regimen.append("对{0}耐药".format(data["var"]["knb"]["evi_sum"]["regimen_R_str"]))
		fjzl_knb_regimen_str = "，".join(fjzl_knb_regimen)

		if data["var"]["knb"] and data["var"]["summary_result"]["format9_forLC10_FJZL"]["drug"]:
			inter_forFJZL.append("本次样本中{0}KRAS/NRAS/BRAF为野生型，提示患者可能{1}。请结合临床。".format(str(data["var"]["summary_result"]["format9_forLC10_FJZL"]["drug"]), str(fjzl_knb_regimen_str)))
		elif data["var"]["knb"]:
			inter_forFJZL.append("本次样本中未检测到靶向用药相关的基因突变，本次样本中KRAS/NRAS/BRAF为野生型，提示患者可能{0}。请结合临床。".format(str(fjzl_knb_regimen_str)))
		elif data["var"]["summary_result"]["format9_forLC10_FJZL"]["drug"]:
			inter_forFJZL.append("本次样本中{0}请结合临床。".format(str(data["var"]["summary_result"]["format9_forLC10_FJZL"]["drug"])))
		else:
			inter_forFJZL.append("本次样本中未检测到靶向用药相关的基因突变。")
		
		# 2025.09.10-新增每个变异结果
		#fjzl_lc10_dataJson = json.dumps({"inter" : inter_forFJZL}, ensure_ascii = False)
		def stran_var_inter(i):
			inter = ""
			if i["bio_category"] == "Snvindel":
				if i["hgvs_p"] != "p.?":
					inter = i["gene_symbol"]+"基因的"+i["hgvs_c"]+" "+i["hgvs_p"]+i["variant_desc_cn"].replace("该变异","")+str(i["variant_interpret_cn"])
				else:
					inter = i["gene_symbol"]+"基因的"+i["hgvs_c"]+i["variant_desc_cn"].replace("该变异","")+str(i["variant_interpret_cn"])
			elif i["bio_category"] == "Cnv":
				inter = "本次实验检测到"+i["gene_symbol"]+"扩增，"+str(i["variant_interpret_cn"])
			elif i["bio_category"] == "Sv":
				if i["variant_interpret_cn"]:
					inter = "本次实验检测到"+i["five_prime_gene"]+"-"+i["three_prime_gene"]+"融合，融合类型为"+i["five_prime_gene"]+":"+i["five_prime_cds"]+\
						"-"+i["three_prime_gene"]+":"+i["three_prime_cds"]+"。"+i["variant_interpret_cn"]
				else:
					inter = "本次实验检测到"+i["five_prime_gene"]+"-"+i["three_prime_gene"]+"融合，融合类型为"+i["five_prime_gene"]+":"+i["five_prime_cds"]+\
						"-"+i["three_prime_gene"]+":"+i["three_prime_cds"]+"。"
			return inter

		def stran_var_for_drug(i):
			inter = ""
			if i["bio_category"] == "Snvindel":
				type_cn = i["type_cn"] if i["type_cn"] != "--" else i["type"]
				type_cn = type_cn.replace("Intronic", "内含子突变").replace("3'UTR","3'UTR突变").replace("5'UTR", "5'UTR突变")
				if i["hgvs_p"] != "p.?":
					inter = "检出"+i["gene_symbol"]+"基因的"+i["hgvs_c"]+" "+i["hgvs_p"]+type_cn
				else:
					inter = "检出"+i["gene_symbol"]+"基因的"+i["hgvs_c"]+type_cn
			elif i["bio_category"] == "Cnv":
				inter = "检出"+i["gene_symbol"]+"扩增"
			elif i["bio_category"] == "Sv":
				inter = "检出"+i["five_prime_gene"]+"-"+i["three_prime_gene"]+"融合"
			return inter

		def fjzl_drug(var):
			var_regimen_S = []
			var_regimen_R = []
			var_regimen = []
			for i in ["A","B","C","D"]:
				if var["var_regimen_forFJZL"]["regimen_"+str(i)+"_S"]:
					var_regimen_S.append("、".join(var["var_regimen_forFJZL"]["regimen_"+str(i)+"_S"])+"（证据等级为"+str(i)+"）")
				if var["var_regimen_forFJZL"]["regimen_"+str(i)+"_R"]:
					var_regimen_R.append("、".join(var["var_regimen_forFJZL"]["regimen_"+str(i)+"_R"])+"（证据等级为"+str(i)+"）")
			if var_regimen_S:
				var_regimen.append("可能对"+"，".join(var_regimen_S)+"敏感")
			if var_regimen_R:
				var_regimen.append("可能对"+"，".join(var_regimen_R)+"耐药")
			return "；".join(var_regimen)
	
		fjzl_var_list = []
		if data["var"]["knb"]:
			knb_regimen = []
			if data["var"]["knb"]["evi_sum"]["regimen_S_str"]:
				knb_regimen.append("对" + data["var"]["knb"]["evi_sum"]["regimen_S_str"]+"敏感")
			if data["var"]["knb"]["evi_sum"]["regimen_R_str"]:
				knb_regimen.append("对" + data["var"]["knb"]["evi_sum"]["regimen_R_str"]+"耐药")
			fjzl_var_list.append(
				{
					"var_name" : "KRAS/NRAS/BRAF WT",
					"interp" : "",
					"popDetail" : "NCCN结直肠癌指南中认为，所有转移性结直肠癌患者均应进行RAS基因（KRAS和NRAS）及BRAF基因突变检测；并指出，KRAS、NRAS基因突变（2/3/4号外显子中）的结直肠癌患者不能用西妥昔单抗或帕尼单抗进行治疗；而BRAF V600E突变的患者除非同时使用BRAF抑制剂，否则极有可能对西妥昔单抗和帕尼单抗无反应。CSCO结直肠癌指南（2022）中，对RAS和BRAF均野生型的患者，推荐使用西妥昔单抗和贝伐珠单抗联合不同化疗；对RAS或BRAF突变的患者，推荐贝伐珠单抗联合不同化疗；对RAS野生、BRAF V600E突变的患者，推荐伊立替康+西妥昔单抗+维莫非尼、BRAF抑制剂+西妥昔单抗±MEK抑制剂。",
					"geneDetail" : "KRAS属于原癌基因，编码的蛋白在细胞内的信号通路中发挥信号转导作用。KRAS突变主要发生在Exon2或3的第12、13和61号密码子上，这些位点位于蛋白的GTP酶结构域。KRAS突变在多种肿瘤中均有发生，包括肺癌，结直肠癌和胰腺癌等等（PMID: 18794081, PMID: 19679400, PMID: 20952405）。NRAS基因属于原癌基因，编码的RAS蛋白在多种细胞信号通路中起着信号转导作用，在细胞的生存与增殖等活动中发挥重要作用。目前，已经在多种肿瘤中发现NRAS基因突变，如黑色素瘤、结直肠癌、甲状腺癌等。BRAF，丝氨酸/苏氨酸蛋白激酶B-raf，为Raf丝氨酸/苏氨酸蛋白激酶家族的成员之一，通过MAP激酶途径发出信号以调节细胞增殖和细胞生长（PMID: 24737949, PMID: 29540830）。已在多种癌症中鉴定出BRAF突变，包括结直肠癌（PMID: 30122982），肺癌（PMID: 29729495），甲状腺（PMID: 12970315）和黑色素瘤（PMID: 24737949），并且还发现了许多突变证明会产生耐药性（PMID: 27478040）。",
					"drugs" : "本次样本中KRAS/NRAS/BRAF为野生型，提示患者可能" + "，".join(knb_regimen)+"。"
				}
			)
		for var in data["var"]["var_somatic"]["level_I"] + \
				   data["var"]["var_somatic"]["level_II"] + \
				   data["var"]["var_somatic"]["level_onco_nodrug"] + \
				   data["var"]["var_somatic"]["level_III"]:
			fjzl_gene_function_list = []
			if var["bio_category"] == "Sv" and "," in var["gene_symbol"] and var["five_prime_gene"] != var["three_prime_gene"]:
				fjzl_gene_function_list.append(var["five_prime_gene_function"])
				fjzl_gene_function_list.append(var["three_prime_gene_function"])
			else:
				fjzl_gene_function_list.append(var["gene_function"])
			fjzl_geneDetail_ = "\n".join(fjzl_gene_function_list)

			popDetail_list = []
			for gene in re.split(",", var["gene_symbol"]):
				if gene in var["var_info_forFJZL"].keys() and var["var_info_forFJZL"][gene]:
					if var["var_info_forFJZL"][gene] not in popDetail_list:
						popDetail_list.append(var["var_info_forFJZL"][gene])
			popDetail_list_ = "\n".join(popDetail_list)

			if var["bio_category"] == "Cnv":
				fjzl_var_list.append(
					{
						"var_name" : var["var_name"] if "var_name" in var.keys() and var["var_name"] else "",
						"interp" : stran_var_inter(var),
						"popDetail" : popDetail_list_,
						"geneDetail" : fjzl_geneDetail_,
						"drugs" : stran_var_for_drug(var) + "，该突变提示" + fjzl_drug(var) + "。" if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"] else "",
						"cn_mean" : var["cn_mean"],
						"cnv_type" : var["cnv_type"]
					}
				)
			else:
				fjzl_var_list.append(
					{
						"var_name" : var["var_name"] if "var_name" in var.keys() and var["var_name"] else "",
						"interp" : stran_var_inter(var),
						"popDetail" : popDetail_list_,
						"geneDetail" : fjzl_geneDetail_,
						"drugs" : stran_var_for_drug(var) + "，该突变提示" + fjzl_drug(var) + "。" if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"] else ""
					}
				)		

		fjzl_data = {}
		fjzl_data["inter"] = inter_forFJZL
		fjzl_data["var"] = fjzl_var_list
		fjzl_lc10_dataJson = json.dumps(fjzl_data, ensure_ascii = False)
		# 2025.09.10-更新完成

		with open(outfile+"/"+json_name+"_fjzl_lc10_inter.json", "w", encoding = "utf-8") as fjzloutFile:
			fjzloutFile.write(fjzl_lc10_dataJson)
	### 福建肿瘤LC10输出结果解读到json结束-2025.07.01

	### 浙江二院小结内容输出到json-2025.10.23
	# 项目包含gBRCA、tBRCA、gHRR、tHRR、ptHRR、150、CP40、MP组织
	# 新增 BPTM Plus组织项目-2026.02.25
	# 新增 BPTM Plus全血项目-2026.04.16
	if data["sample"]["company"] == "浙江大学医学院附属第二医院" and (re.search("ZJFE", report_name) or re.search("ZJEY", report_name)) and data["sample"]["prod_names"] in \
	["BRCA1/BRCA2（全血）", "BRCA1/BRCA2（组织）", "HRR（全血）", "HRR（组织）", "HRR（组织 全血）", "遗传易感150基因", "Classic Panel", "Master Panel（组织）", "BPTM Plus（组织）", "BPTM Plus（全血）"]:
		zjfe_sum_data = get_summary_for_specialcomany.ZJFE_summary(data, judge_brca_cnv)
		sum_dataJson = json.dumps({"summary" : zjfe_sum_data}, ensure_ascii = False)
		with open(outfile+"/"+json_name+"_zjfe_summary.json", "w", encoding = "utf-8") as outFile:
			outFile.write(sum_dataJson)
	### 浙江二院小结内容输出到json结束-2025.10.23

	# 参数为T时，填充用数据转化为json输出，便于开发
	if outjson == "T":
		dataJson = json.dumps(data, ensure_ascii = False)
		with open(outfile+"/"+json_name+"_to_word.json", "w", encoding = "utf-8") as outFile:
			outFile.write(dataJson)


	# 加一个分页符-2023.05.17
	data["page_break"] = R("\f")
	# 加个判断字段，为mlpa时报告展示MLPA，为gcnv时展示gCNV-2024.08.30
	#data["judge_brca_cnv"] = "mlpa"
	# 模板填充
	print (report_name)
	print (judge_brca_cnv)
	if report_name:
		path = os.path.join(report_template, "template_main", report_name)
		tpl = DocxTemplate(path)
		# 图片填充
		# 2025.09.29-多产品出一份报告的图片在定制代码里加，这边为通用的
		if not ("merge_order" in jsonDict.keys() and jsonDict["merge_order"] and len(jsonDict["merge_order"]["cnv"]["sample_id_list"]) >= 2 and report_name_judge_merge):
			data["image"] = getImage.render_image(tpl, data, jsonDict, report_name, image, config)
		print ("填充数据生成完毕！", datetime.datetime.now())
		print ("开始填充报告：", datetime.datetime.now())
		# 拼接模板-20221216
		# 2025.09.29-基础模板不兼容多产品出一份报告，这边为通用的
		if not ("merge_order" in jsonDict.keys() and jsonDict["merge_order"] and len(jsonDict["merge_order"]["cnv"]["sample_id_list"]) >= 2 and report_name_judge_merge):
			data["subdoc"] = baseTemplate.BaseReport(data, tpl, merge_template, report_template, json_name, outfile)
		tpl.render(data)
		tpl.save(outfile+"/"+json_name+".docx")
		print ("报告填充完成！", datetime.datetime.now())
		

		# 加个自动更新域功能-仅供需要更新目录的模板使用，待开发-20220922
		#judge_update = "yes"
		#if judge_update:
			#namespace = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
			#element_updatefields = lxml.etree.SubElement(tpl.settings.element, namespace+"updateFields")
			#element_updatefields.set(namespace+"val", "true")
			#element_updatefields.set(namespace+"val","true")
			#tpl.save(outfile+"/"+json_name+".docx")
	else:
		print ("未匹配到报告模板")
	



def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("-s", "--json_name", dest = "json_name", required = True)
	parser.add_argument("-o", "--outfile", dest = "outfile", required = True)
	parser.add_argument("-c", "--config", dest = "config", required = True)
	parser.add_argument("-r", "--report_template", dest = "report_template", required = True)
	parser.add_argument("-j", "--outjson", dest = "outjson", required = True)
	parser.add_argument("-i", "--image", dest = "image", required = True)
	arg = parser.parse_args()

	return arg

if __name__ == '__main__':
	args = parse_args()
	get_data(json_name=args.json_name, outfile=args.outfile, config=args.config, report_template=args.report_template, outjson=args.outjson, image = args.image)