#-*- coding:gbk -*-
import datetime
import copy
import re

'''
Discription
	
	获取json文件中的sample info，转化为报告模板方便填充的格式。 

'''

def getSample(jsonDict):
	data = copy.deepcopy(jsonDict["sample_info"])
	data["report_date"] = str(datetime.date.today())
	for k, v in data.items():
		if not v:
			data[k] = ""
	# 2026.01.30-接收日期更新
	# LIMS临检订单，接收日期分为组织样本接收日期tissue_date_received和血液样本接收日期blood_date_received 
	# 单样本产品和血液配对产品，组织项目把tissue_date_received传递给receive_data，血液项目把blood_date_received传递给receive_data（报告模板改动最小）
	# 组织-血液配对样本，则把tissue_date_received传递给receive_data，对照样本的接收日期在模板里更改。
	if "order_type" in data.keys() and data["order_type"] and data["report_module_type"] == "rummage":
		single_blood_prod = ["10基因（血液）", "61遗传基因", "HRR（全血）", "Pan116（血液）", "BRCA1/BRCA2（全血）", \
					   		 "林奇综合征", "TC21（血液）", "GA18（血液）", "LC76（血液）", "CRC25（血液）", \
							 "BPTM（全血）", "BRCA1/2（扩增子）", "Master Panel（血液）", "遗传易感150基因", "Master Panel（血液-综合）", \
							 "Pan116（血液）-胰腺", "Pan116（血液）-乳腺", "HRR（全血）-前列腺", "BPTM Plus（全血）", "Pan116（血液）-黑色素瘤", \
							 "Pan116（血液）-头颈", "PTM Plus（全血）-林奇", "OncoPro（血液）", "遗传易感150基因-乳腺癌-华西", "遗传易感150基因-华西-体检", \
							 "OncoPro Panel（血液）-5基因", "遗传易感150基因-前列腺癌-华西", "乳腺癌血浆5基因", "遗传易感150基因-华西-消化道肿瘤", "OncoPro Panel（血液）", "MRD（LC10）"]
		if data["prod_names"] in single_blood_prod:
			data["receive_data"] = data["blood_date_received"] if "blood_date_received" in data.keys() and data["blood_date_received"] else data["receive_data"]
			print ("LIMS临检订单，血液样本接收日期", data["receive_data"])
		else:
			data["receive_data"] = data["tissue_date_received"] if "tissue_date_received" in data.keys() and data["tissue_date_received"] else data["receive_data"]
			print ("LIMS临检订单，组织样本接收日期", data["receive_data"])
	# 2026.01.30-接收日期更新结束

	# 2025.08.20-富集前/后肿瘤细胞含量字段，v4没有数据的时候会返回“N/A”，这个兼容下
	data["tumor_content"] = data["tumor_content"] if "tumor_content" in data.keys() and data["tumor_content"] and data["tumor_content"] != "N/A" else ""
	data["tumor_content_macrodissection_performed"] = data["tumor_content_macrodissection_performed"] if "tumor_content_macrodissection_performed" in data.keys() and data["tumor_content_macrodissection_performed"] and data["tumor_content_macrodissection_performed"] != "N/A" else ""
	# 2025.08.20
		
	# 2025-08-15-肿瘤细胞含量返回规则有更新
	# v3-富集前/富集后由报告系统判定，返回数值大的那个到tumor_content
	# v4-富集前tumor_content/富集后tumor_content_macrodissection_performed分开字段返回
	# 体系和胚系肿瘤细胞含量展示规则不同
	# 体系：有富集后肿瘤细胞含量就展示富集后的，没有就展示富集前的
	# 胚系：报告中富集前/后结果都要展示，但如果有富集后结果，报告里阈值判定要看富集后，否则就看富集前（模板中也要调整）
	somatic_prod_list = ["10基因（组织）", "CRC12-MSI", "Pan116（组织）", "Master Panel（组织）", "Classic Panel", "TC21（组织）", \
					  	 "GA18（组织）", "LC76（组织）", "CRC25（组织）", "Master Panel（组织-综合）", "Classic Panel-胃", "Classic Panel-胃肠", \
						 "Classic Panel-甲状腺", "Classic Panel-肝胆", "Classic Panel-骨", "Classic Panel-黑色素", "Classic Panel-头颈", \
						 "Pan116（组织）-胰腺", "Pan116（组织）-乳腺", "OncoPro（组织）", "Master Panel V2（组织）", "Master Panel V1（组织）", \
						 "RNASeq（肉瘤）", "Pan116（组织）-黑色素瘤", "Pan116（组织）-头颈", "OncoPro（组织）-妇科肿瘤", "Master Panel V1（组织_北医）", \
						 "Classic Panel 200（组织）"]
	if data["prod_names"] in somatic_prod_list:
		data["tumor_content"] = data["tumor_content_macrodissection_performed"] if "tumor_content_macrodissection_performed" in data.keys() and data["tumor_content_macrodissection_performed"] else data["tumor_content"]
	# 2025-08-15-新增完成
	
	#新增日期格式：XX年XX月XX日
	receive_date = re.split("-", data["receive_data"])
	data["receive_date_special_1"] = receive_date[0] + "年" + receive_date[1] + "月" + receive_date[2] + "日" if receive_date and len(receive_date) >= 3 else ""
	report_date = re.split("-", data["report_date"])
	data["report_date_special_1"] = report_date[0] + "年" + report_date[1] + "月" + report_date[2] + "日" if report_date and len(report_date) >= 3 else ""

	#新增生信分析时间
	json_name_list = re.split("_", data["json_batch_name"])
	json_date = json_name_list[0] if json_name_list else ""
	data["json_date"] = json_date[:4] + "年" + json_date[4:6] + "月" + json_date[6:] + "日" if json_date and len(json_date) == 8 else ""

	#判断临检样本是厦门接收还是上海接收的，以首字母是否带S进行判断
	data["locate"] = "SH" if re.match("S", data["sample_id"]) else "XM"

	# mark和note字段改了2025.05.08
	# BIMS传过来的寄送备注输出到note中，配置输出json为mark；
	# 手动上传以后改用备注note，配置输出为mark
	# 原来的寄送备注mark就不用了，输出为note
	mark = data["mark"] if "mark" in data.keys() and data["mark"] else ""
	note = data["note"] if "note" in data.keys() and data["note"] else ""
	mark_and_note_list = []
	if mark:
		mark_and_note_list.append(mark)
	if note:
		mark_and_note_list.append(note)
	data["mark"] = "；".join(mark_and_note_list)

	# 浙江妇保，寄送备注处理
	data["ZJFB_mark_dict"] = {}
	if data["mark"]:
		for i in re.split("；|;", data["mark"]):
			tmp = re.split("：|:", i)
			if len(tmp) == 2:
				data["ZJFB_mark_dict"][tmp[0]] = ""
				data["ZJFB_mark_dict"][tmp[0]] = tmp[1]
	
	# 浙江妇保，寄送备注处理
	data["ZJFB_mark_dict"] = {}
	if data["mark"]:
		for i in re.split("；|;", data["mark"]):
			tmp = re.split("：|:", i)
			if tmp:
				data["ZJFB_mark_dict"][tmp[0]] = tmp[1] if len(tmp) >= 2 else ""


			if len(tmp) == 2:
				data["ZJFB_mark_dict"][tmp[0]] = ""
				data["ZJFB_mark_dict"][tmp[0]] = tmp[1]
	
	# 肿瘤细胞含量新增数值字段，便于比对
	data["tumor_content_num"] = data["tumor_content"].replace("%", "") if "tumor_content" in data.keys() and data["tumor_content"] else ""

	# 2025-08-15富集后肿瘤细胞含量新增数值字段，便于比对
	data["tumor_content_macrodissection_performed_num"] = data["tumor_content_macrodissection_performed"].replace("%", "") if "tumor_content_macrodissection_performed" in data.keys() and data["tumor_content_macrodissection_performed"] else ""

	# 特殊处理：临检项目，年龄为0的放空，性别为未知的放空-20220923
	if jsonDict["sample_info"]["report_module_type"] != "hospital":
		data["age"] = "" if data["age"] and data["age"] == "0" else data["age"]
		data["gender"] = "" if data["gender"] and data["gender"] == "未知" else data["gender"]

	# 香港大学深圳医院，寄送备注处理-2022.10.10
	data["XGDX_mark_dict"] = {}
	if data["mark"]:
		for i in re.split("，", data["mark"]):
			tmp = re.split("：", i)
			if tmp:
				data["XGDX_mark_dict"][tmp[0]] = tmp[1] if len(tmp) >= 2 else ""
	
	# 采集日期，组织采集日期还是原来的字段：gather_data
	# 新增一个血液采集日期，这边做兼容-2022.11.29
#	data["blood_collection_date"] = data["blood_collection_date"] if "blood_collection_date" in data.keys() and data["blood_collection_date"] else ""
	data["blood_collection_date"] = re.split(" ", data["blood_collection_date"])[0] if "blood_collection_date" in data.keys() and data["blood_collection_date"] else ""

	# 日期删除具体时间这部分代码注释掉，有些医院要展示具体时间(北京协和BRCA)，有需要删除具体时间的话，由报告系统处理吧-2023.12.19
	# 现在日期返回格式有的会带具体时间，这边加个兼容处理-刘炜芬-20231025
	# 2024.03.06-增加blood_date_received和tissue_date_received（药企模板使用）
	date_key = ["receive_data", "tissue_collection_date", "blood_collection_date", "gather_data", "submission_date", \
				"application_date", "sampling_time", "analysis_date", "section_date", "blood_date_received", "tissue_date_received"]
	for key in date_key:
		if key in data.keys():
			data[key] = re.split(" ", data[key])[0] if data[key] else ""
	# 添加完成-20231025

	# 新增兼容-组织采集日期，目前大部分模板都是使用	gather_data字段，新增了一个tissue_collection_date-刘炜芬，2023.11.14
	# 若存在gather_data，则gather_data使用该字段内容，若不存在gather_data，则使用tissue_collection_date
	data["gather_data"] = data["gather_data"] if "gather_data" in data.keys() and data["gather_data"] else \
						  data["tissue_collection_date"] if "tissue_collection_date" in data.keys() and data["tissue_collection_date"] else \
						  ""
	# 兼容完成-2023.11.14

	# 新增需求-2024.01.04
	# 北京大学人民医院-临检，送蜡卷的话没有肿瘤细胞含量，销售写在寄送备注中，由报告脚本提取；
	# 寄送备注格式：肿瘤细胞含量：XX%；XXXXX。
	# 填写“XX%”时若模板中需要综合评估，可直接判断；写为“>30%”的，需要模板中做兼容，兼容的话考虑定制的就行，临检通用模板都没有综合评估了。
	if data["company"] == "北京大学人民医院" and data["report_module_type"] == "rummage":
		bdrm_mark = {}
		if data["mark"]:
			for i in re.split(";|；", data["mark"]):
				tmp = re.split("：|:", i)
				if tmp:
					bdrm_mark[tmp[0]] = tmp[1] if len(tmp) >= 2 else ""
		if not data["tumor_content"]:
			data["tumor_content"] = bdrm_mark.get("肿瘤细胞含量", "")
		#print (bdrm_mark)
	# 更新完成-2024.01.04

	# 2025.09.08-增加一个检测平台，临检通用模板中的检测平台统一从这个字段中获取，方便后续更新
	#data["ngs_platform"] = "Illumina/贝瑞基因测序平台、ADx-SEQ200 Plus测序仪"
	# 增加赛陆测序仪，2025.11.20正式启用
	data["ngs_platform"] = "Illumina/贝瑞基因测序平台、ADx-SEQ200 Plus测序仪、赛陆测序仪"
	# 2025.09.08-增加完成

	# 2025.12.09-CP43部分进院模板启用通用最新版-这边加个启用医院列表
	# 2026.03.18-增加华东医院
	# 2026.03.24-增加衡水市第二人民医院
	# 2026.04.01-增加上海交通大学医学院附属新华医院
	# 2026.04.02-衡水市第二人民医院改为40基因
	# 2026.04.07-吉林大学第二医院改为43基因
	# 2026.04.10-聊城市人民医院改为43基因
	# 2026.04.16-衡水市第二人民医院改为43基因
	data["cp43_special_company_list"] = ["四川大学华西医院", "德阳市人民医院", "华东医院", "上海交通大学医学院附属新华医院", "吉林大学第二医院", "聊城市人民医院", "衡水市第二人民医院"]
	# 2025.12.09-增加完成

	# 2026.02.27-协会名称写在这，模板中调用，有修改的话就不用逐份更新了
	data["association"] = {}
	data["association"]["AMP"] = "美国分子病理学协会（Association for Molecular Pathology，AMP）"
	data["association"]["ASCO"] = "美国临床肿瘤学会（American Society of Clinical Oncology，ASCO）"
	data["association"]["CAP"] = "美国病理学家协会（College of American Pathologists，CAP）"
	data["association"]["ACMG"] = "美国医学遗传学和基因组学学会（American College of Medical Genetics and Genomics，ACMG）"
	data["association"]["ClinGen"] = "临床基因组资源中心（Clinical Genome Resource，ClinGen）"
	data["association"]["CGC"] = "癌症基因组学联盟（Cancer Genomics Consortium，CGC）"
	data["association"]["VICC"] = "癌症变异解读联盟（Variant Interpretation for Cancer Consortium，VICC）"
	data["association"]["IARC"] = "国际癌症研究所（International Agency of Research on Cancer，IARC）"
	data["association"]["AMP_cn"] = "美国分子病理学协会"
	data["association"]["ASCO_cn"] = "美国临床肿瘤学会"
	data["association"]["CAP_cn"] = "美国病理学家协会"
	data["association"]["ACMG_cn"] = "美国医学遗传学和基因组学学会"
	data["association"]["ClinGen_cn"] = "临床基因组资源中心"
	data["association"]["CGC_cn"] = "癌症基因组学联盟"
	data["association"]["VICC_cn"] = "癌症变异解读联盟"
	data["association"]["IARC_cn"] = "国际癌症研究所"
	# 2026.02.28-部分需要带版本的指南也写在这，模板调用
	data["guideline"] = {}
	data["guideline"]["germline_var"] = "《遗传变异分类标准与指南》（2015年版）"
	data["guideline"]["somatic_var"] = "《癌症中体细胞变异致病性分类标准（致癌性）》（2022年版）"
	data["guideline"]["tumor_var"] = "《肿瘤变异解读及报告指南》（2017年版）"
	data["guideline"]["ngs_clinic_report_inter"] = "《二代测序临床报告解读指引》"
	# 2026.02.27-新增完成
	
	return data