#-*- coding:gbk -*-
#import os
#import xlrd
import re
from libs.getConfig import getconfigxlsx
from libs.rule import judge_var, judgeRegimen
from libs.getConfig import get_io_refer
from libs.getConfig import get_io_zsly
from customize_filters import fdzs_mp_cnv_stran

'''
Discription
	
	免疫正负相关检测结果。 
	适用：
	1. v3临检通用MP IVD
	2. v3临检通用116
	3. v3浙肿Master 跟临检的格式不一样，为了防止出错，额外再处理一下
	4. 新增北京医院-20231120
	5. 新增CP200-2024.10.10
	
'''
def get_io_detect(var_data, config):
	io = {}
	io["result"], io["io_p_summary"], io["io_n_summary"] = io_detect(var_data)
	io["result_116"], io["io_p_summary_116"], io["io_n_summary_116"] = io_detect_for_116(var_data)
	io["io_p_summary_ZJZL"], io["io_n_summary_ZJZL"] = io_detect_for_ZJZL(var_data, config)
	io["result_new_cnvlist"], io["io_p_summary_new_cnvlist"], io["io_n_summary_new_cnvlist"] = io_detect_new_cnvlist(var_data)
	# 新增北京医院-合并梦晨代码-2023.11.20
	io['io_BJYY'], io['io_p_BJYY'], io['io_n_BJYY'], io['io_BJYY_num'] = io_detect_for_BJYY(var_data)
	# 合并结束-2023.11.20
	# 2025.07.21-北京医院新增HD
	io['io_BJYY_hd'], io['io_p_BJYY_hd'], io['io_n_BJYY_hd'], io['io_BJYY_num_hd'] = io_detect_for_BJYY_hd(var_data)
	# 2025.07.21-新增完成
	# 新增孙逸仙116-三字母-2024.08.19
	io["result_116_syx"], io["io_p_summary_116_syx"], io["io_n_summary_116_syx"] = io_detect_for_116_syx(var_data)
	# 新增结束-2024.08.19
	# 新增CP200-2024.10.10
	io["io_p_summary_cp200"], io["io_n_summary_cp200"] = io_detect_for_cp200(var_data)
	# 新增中山六院，删除CDKN2B基因-2025.04.09
	io["io_p_summary_cp200_zsly"], io["io_n_summary_cp200_zsly"] = io_detect_for_cp200_zsly(var_data)
	# 新增结束-2024.10.10
	# 新增oncopro血液-2024.12.19
	# 新增返回io_result-2026.01.12
	#io["io_p_summary_boncopro"], io["io_n_summary_boncopro"] = io_detect_for_boncopro(var_data)
	io["result_boncopro"], io["io_p_summary_boncopro"], io["io_n_summary_boncopro"] = io_detect_for_boncopro(var_data)

	# 2025.02.21-返回免疫正负相关参考文献-从配置表中获取
	io["io_refer"] = get_io_refer(config)

	# 2025.03.24-返回免疫正负相关-癌种-解析-从配置表中获取
	io["tumor_inter"] = get_io_zsly(config)

	# 2025.04.02-中山六院-io结果需要带频率等信息
	io["zsly_result"] = ZSLY_io_detect(var_data)

	# 2025.06.17-新增展示HD的结果
	io["hd_result"], io["hd_io_p_summary"], io["hd_io_n_summary"] = io_detect_hd(var_data)

	# 2025.07.17-新增展示HD结果
	io["io_p_summary_ZJZL_hd"], io["io_n_summary_ZJZL_hd"] = io_detect_for_ZJZL_hd(var_data, config)

	# 2026.02.03-新增武汉协和-氨基酸展示三字母
	io["abbr_hd_result"], io["abbr_hd_io_p_summary"], io["abbr_hd_io_n_summary"] = abbr_io_detect_hd(var_data)

	# 2026.02.25-新增北三MP-胚系3类有用药位点不展示
	io["hd_result_bds"], io["hd_io_p_summary_bds"], io["hd_io_n_summary_bds"] = io_detect_hd_bds(var_data)

	# 2026.03.17-新增北大人民MP，格式与北京医院不同
	io['io_bdrm'], io['io_p_bdrm'], io['io_n_bdrm'], io['io_bdrm_num'] = io_detect_for_bdrm(var_data)
	# 2026.05.08-新增北大人民MP，包含HD结果
	io['io_bdrm_hd'], io['io_p_bdrm_hd'], io['io_n_bdrm_hd'], io['io_bdrm_num_hd'] = io_detect_for_bdrm_hd(var_data)
	# 2026.05.14 - 复旦中山MP，相关变异CNV描述修改，MET/ERBB2 5-9报拷贝数增加，≥10报扩增；其他基因5-14报拷贝数增加，≥15报扩增，嵇梦晨，2026.05.14
	io['io_fdzs_mp'], io['io_p_fdzs_mp'], io['io_n_fdzs_mp'] = io_detect_fdzs_mp(var_data)

	return io

def io_detect(var_data):
	# 返回结果中的io_result用于填充IO表，", ".join(io_p_list), ", ".join(io_n_list)用于填充检测结果小结
	io_result = {}
	# 汇总体细胞I/II/肿瘤发生发展相关变异+胚系致病/疑似致病性变异
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","ERCC1","FANCA","MRE11",\
				 "PALB2","RAD50","XRCC1","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53",\
				 "KRAS","CD274","ARID1A","LRP1B","SETD2","PRKDC","TERT","KMT2D","FAT1","CDK12",\
				 "SERPINB3","SERPINB4"]
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1","IRF1",\
				 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","CCND1","FGF3","FGF4","FGF19"]
	cnv_gene_list = ["CD274", "MDM2", "MDM4", "CCND1", "FGF3", "FGF4", "FGF19"]

	level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]
	# 2025.02.21- 更新展示变异等级
	# 1. 有配对：体细胞I/II/肿瘤发生发展 + 胚系4/5 + 胚系3类但有用药
	# 2. 无配对：体细胞I/II/肿瘤发生发展 + 预测胚系4/5（归为I/II/肿瘤发生发展） + 预测胚系3类但有用药
	# ==> 即新增（确认/预测）胚系3类但有用药的变异即可
	germline_level3_regimen = [var for var in var_data if var["var_origin"] == "germline" and var["clinic_num_g"] == 3 and judgeRegimen(var)]
	if germline_level3_regimen:
		level_12_var.extend(germline_level3_regimen)
	# 2025.02.21-新增完成


	for var in level_12_var:
		# 仅展示扩增的基因
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in cnv_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append("扩增")
		# 仅展示融合的基因
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			# 融合可能会出现exon相同断点不同的变异，报告中会重复，这边加个去重-2023.04.13
			if var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"融合" not in io_result["ALK"]:
				io_result["ALK"].append(var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"融合")
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["hgvs_p"] != "p.?":
				io_result[var["gene_symbol"]].append(var["hgvs_p"])
			else:
				io_result[var["gene_symbol"]].append(var["hgvs_c"])
	# summary展示
	both_cnv_list = ["CCND1","FGF3","FGF4","FGF19"]
	io_p_list = ["{0} {1}".format(k, i) for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_P]
	io_n_list = ["{0} {1}".format(k, i) if not re.search("融合", i) else i for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_N]
	#io_p_list = []
	#io_n_list = []
	#for k, v in io_result.items():
	#	io_p_list.extend(["{0} {1}".format(k, i) for i in v if k not in both_cnv_list and k in io_gene_P])
	#	io_n_list.extend(["{0} {1}".format(k, i) for i in v if k not in both_cnv_list and k in io_gene_N])


	# 处理CNV共突变
	if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF4" in io_result.keys() and "FGF19" in io_result.keys():
		io_n_list.append("CCND1/FGF3/FGF4/FGF19扩增")
		#print (", ".join(io_n_list))
	return io_result, ", ".join(io_p_list), ", ".join(io_n_list)

def io_detect_new_cnvlist(var_data):
	# 返回结果中的io_result用于填充IO表，", ".join(io_p_list), ", ".join(io_n_list)用于填充检测结果小结
	io_result = {}
	# 汇总体细胞I/II/肿瘤发生发展相关变异+胚系致病/疑似致病性变异
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","ERCC1","FANCA","MRE11",\
				 "PALB2","RAD50","XRCC1","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53",\
				 "KRAS","CD274","ARID1A","LRP1B","SETD2","PRKDC","TERT","KMT2D","FAT1","CDK12",\
				 "SERPINB3","SERPINB4"]
	#io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1","IRF1",\
	#			 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","CCND1","FGF3","FGF4","FGF19"]
	#cnv_gene_list = ["CD274", "MDM2", "MDM4", "CCND1", "FGF3", "FGF4", "FGF19"]
	# 删除FGF4扩增-2022307.20
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1","IRF1",\
				 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","CCND1","FGF3","FGF19"]
	cnv_gene_list = ["CD274", "MDM2", "MDM4", "CCND1", "FGF3", "FGF19"]
	
	level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]
	# 2025.02.21- 更新展示变异等级
	# 1. 有配对：体细胞I/II/肿瘤发生发展 + 胚系4/5 + 胚系3类但有用药
	# 2. 无配对：体细胞I/II/肿瘤发生发展 + 预测胚系4/5（归为I/II/肿瘤发生发展） + 预测胚系3类但有用药
	# ==> 即新增（确认/预测）胚系3类但有用药的变异即可
	germline_level3_regimen = [var for var in var_data if var["var_origin"] == "germline" and var["clinic_num_g"] == 3 and judgeRegimen(var)]
	if germline_level3_regimen:
		level_12_var.extend(germline_level3_regimen)
	# 2025.02.21-新增完成

	for var in level_12_var:
		# 仅展示扩增的基因
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in cnv_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append("扩增")
		# 仅展示融合的基因
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			# 融合可能会出现exon相同断点不同的变异，报告中会重复，这边加个去重-2023.04.13
			if var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"融合" not in io_result["ALK"]:
				io_result["ALK"].append(var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"融合")
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["hgvs_p"] != "p.?":
				io_result[var["gene_symbol"]].append(var["hgvs_p"])
			else:
				io_result[var["gene_symbol"]].append(var["hgvs_c"])
	# summary展示
	#both_cnv_list = ["CCND1","FGF3","FGF4","FGF19"]
	# 删除FGF4扩增-2022307.20
	both_cnv_list = ["CCND1","FGF3","FGF19"]
	io_p_list = ["{0} {1}".format(k, i) for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_P]
	io_n_list = ["{0} {1}".format(k, i) if not re.search("融合", i) else i for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_N]
	#io_p_list = []
	#io_n_list = []
	#for k, v in io_result.items():
	#	io_p_list.extend(["{0} {1}".format(k, i) for i in v if k not in both_cnv_list and k in io_gene_P])
	#	io_n_list.extend(["{0} {1}".format(k, i) for i in v if k not in both_cnv_list and k in io_gene_N])


	# 处理CNV共突变
	#if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF4" in io_result.keys() and "FGF19" in io_result.keys():
	#	io_n_list.append("CCND1/FGF3/FGF4/FGF19扩增")
	# 删除FGF4扩增-2022307.20
	if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF19" in io_result.keys():
		io_n_list.append("CCND1/FGF3/FGF19扩增")
	#print (", ".join(io_n_list))
	return io_result, ", ".join(io_p_list), ", ".join(io_n_list)

def io_detect_for_116(var_data):
	# 返回结果中的io_result用于填充IO表，", ".join(io_p_list), ", ".join(io_n_list)用于填充检测结果小结
	io_result = {}
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","FANCA","MRE11","PALB2","MLH1","MSH2","MSH6",\
				 "PMS2","POLE","POLD1","TP53","KRAS","CD274","ARID1A","TERT","CDK12"]
	#io_gene_N = ["EGFR","ALK","CDKN2A","CDKN2B","STK11","JAK1","JAK2","APC","CTNNB1","PTEN"]
	# 新增CCND1/FGF3/FGF19共扩增-2023.07.12
	io_gene_N = ["EGFR","ALK","CDKN2A","CDKN2B","STK11","JAK1","JAK2","APC","CTNNB1","PTEN", "CCND1", "FGF3", "FGF19"]
	level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]
	# 2025.02.21- 更新展示变异等级
	# 1. 有配对：体细胞I/II/肿瘤发生发展 + 胚系4/5 + 胚系3类但有用药
	# 2. 无配对：体细胞I/II/肿瘤发生发展 + 预测胚系4/5（归为I/II/肿瘤发生发展） + 预测胚系3类但有用药
	# ==> 即新增（确认/预测）胚系3类但有用药的变异即可
	germline_level3_regimen = [var for var in var_data if var["var_origin"] == "germline" and var["clinic_num_g"] == 3 and judgeRegimen(var)]
	if germline_level3_regimen:
		level_12_var.extend(germline_level3_regimen)
	# 2025.02.21-新增完成
	for var in level_12_var:
		# 仅展示扩增的基因
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in ["CD274", "CCND1", "FGF3", "FGF19"]:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append("扩增")
		# 仅展示融合的基因
		if var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			io_result["ALK"].append(var["five_prime_gene"]+"-"+var["three_prime_gene"]+"融合")
			# 加一个孙逸仙的 Gene1-ALK(G8:A12)基因重排 -2023.08.02
			if "ALK_syx" not in io_result.keys():
				io_result.setdefault("ALK_syx", [])
			io_result["ALK_syx"].append("{0}-{1}({2}{3}:{4}{5})基因重排".format(var["five_prime_gene"], var["three_prime_gene"], var["five_prime_gene"][0], \
								var["five_prime_cds"].replace("exon", ""), var["three_prime_gene"][0], var["three_prime_cds"].replace("exon", "")))
			# 新增结束-2023.08.02
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["hgvs_p"] != "p.?":
				io_result[var["gene_symbol"]].append(var["hgvs_p"])
			else:
				io_result[var["gene_symbol"]].append(var["hgvs_c"])
	# summary展示
	both_cnv_list = ["CCND1", "FGF3", "FGF19"]
	io_p_list = ["{0} {1}".format(k, i) for k,v in io_result.items() for i in v if k in io_gene_P]
	io_n_list = ["{0} {1}".format(k, i) if not re.search("融合", i) else i for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_N]
	# 处理CNV共突变
	if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF19" in io_result.keys():
		io_n_list.append("CCND1/FGF3/FGF19扩增")
	
	return io_result, ", ".join(io_p_list), ", ".join(io_n_list)

def io_detect_for_ZJZL(var_data, config):
	fdzs_dict, fjzl_database, Data = getconfigxlsx(config)
	# io_result用于填充IO表
	io_result = {}

	# 删除FGF4扩增-2025.02.21
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","ERCC1","FANCA","MRE11",\
				 "PALB2","RAD50","XRCC1","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53",\
				 "KRAS","CD274","ARID1A","LRP1B","SETD2","PRKDC","TERT","KMT2D","FAT1","CDK12",\
				 "SERPINB3","SERPINB4"]
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1","IRF1",\
				 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","CCND1","FGF3","FGF19"]
	cnv_gene_list = ["CD274", "MDM2", "MDM4", "CCND1", "FGF3", "FGF19"]
	# 展示I/II类和肿瘤发生发展相关变异
	level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]
	# 2025.02.21- 更新展示变异等级
	# 1. 有配对：体细胞I/II/肿瘤发生发展 + 胚系4/5 + 胚系3类但有用药
	# 2. 无配对：体细胞I/II/肿瘤发生发展 + 预测胚系4/5（归为I/II/肿瘤发生发展） + 预测胚系3类但有用药
	# ==> 即新增（确认/预测）胚系3类但有用药的变异即可
	germline_level3_regimen = [var for var in var_data if var["var_origin"] == "germline" and var["clinic_num_g"] == 3 and judgeRegimen(var)]
	if germline_level3_regimen:
		level_12_var.extend(germline_level3_regimen)
	# 2025.02.21-新增完成
	for var in level_12_var:
		# 仅展示扩增的基因
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in cnv_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append(
				{
					"var_type" : "cnv", 
					"var_info" : "扩增"
					}
				)
		# 仅展示融合的基因
		# 可能存在两个融合基因都在检测范围里的情况，系统json返回的gene_symbol为gene1,gene2，需额外做识别
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			io_result["ALK"].append(
				{
					"var_type" : "sv", 
					"three_prime_gene" : var["three_prime_gene"], 
					"three_prime_cds" : var["three_prime_cds"], 
					"five_prime_gene" : var["five_prime_gene"], 
					"five_prime_cds" : var["five_prime_cds"]
					}
				)
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append(
				{
					"var_type" : "snvindel", 
					"hgvs_p" : var["hgvs_p"], 
					"hgvs_c" : var["hgvs_c"], 
					"var_origin" : var["var_origin"], 
					"gene_region" : var["gene_region"], 
					"transcript_primary" : var["transcript_primary"]
					}
				)
	# summary展示
	io_inter = []
	for k, v in io_result.items():
		if k not in ["CCND1","FGF3","FGF19"]:
			if k in Data.keys():
				io_inter.append({
					"gene_symbol" : k,
					"var_info" : v,
					"inter" : Data[k]
				})
	if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF19" in io_result.keys():
		io_inter.append({
			"gene_symbol" : "CCND1/FGF3/FGF19",
			"var_info" : io_result["CCND1"],
			"inter" : Data["CCND1"] 
		})
	
	io_p_list = [i for i in io_inter if i["gene_symbol"] in io_gene_P]
	io_n_list = [i for i in io_inter if i["gene_symbol"] in io_gene_N or i["gene_symbol"] == "CCND1/FGF3/FGF19"]

	return io_p_list, io_n_list

# 北京医院-合并梦晨代码-2023.11.20
def io_detect_for_BJYY(var_data):
	# io_result用于填充IO表
	io_result = {}
	# io_p_result用于填充检测结果小结
	io_p_result = ""
	# io_n_result用于填充检测结果小结
	io_n_result = ""
	# 汇总体细胞I/II/肿瘤发生发展相关变异+胚系致病/疑似致病性变异
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","ERCC1","FANCA","MRE11","PALB2","RAD50","XRCC1","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53","KRAS","CD274","ARID1A","LRP1B","SETD2","PRKDC","TERT","KMT2D","FAT1","CDK12","SERPINB3","SERPINB4"]
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1","IRF1","JAK1","JAK2","APC","CTNNB1","B2M","PTEN","CCND1","FGF3","FGF19"]
	level_12_var = [var for var in var_data if (var["var_origin"] != "germline" and var["clinic_num_s"] in [5, 4]) or (var["var_origin"] == "germline" and var["clinic_num_g"] in [5, 4])]
	# 2025.02.21- 更新展示变异等级
	# 1. 有配对：体细胞I/II/肿瘤发生发展 + 胚系4/5 + 胚系3类但有用药
	# 2. 无配对：体细胞I/II/III类 + 预测胚系4/5（归为肿瘤发生发展） + 预测胚系3类但有用药
	# ==> 即新增（确认/预测）胚系3类但有用药的变异即可
	germline_level3_regimen = [var for var in var_data if var["var_origin"] == "germline" and var["clinic_num_g"] == 3 and judgeRegimen(var)]
	if germline_level3_regimen:
		level_12_var.extend(germline_level3_regimen)
	# 2025.02.21-新增完成

	for var in level_12_var:
		# 仅展示扩增的基因
		# 扩增共突变
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in ["CD274", "MDM2", "MDM4", "CCND1", "FGF3", "FGF19"]:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append("扩增")
		# 仅展示融合的基因
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			io_result["ALK"].append(var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"融合")
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["hgvs_p"] != "p.?":
				io_result[var["gene_symbol"]].append(var["hgvs_p_BJYY"])
			else:
				io_result[var["gene_symbol"]].append(var["hgvs_c"])
	
	# summary展示
	io_p_list_BJYY = []
	io_n_list_BJYY = []
	for k, v in io_result.items():
		if k not in ["CCND1","FGF3","FGF19"]:
			for i in v:
				if k in io_gene_P:
					io_p_list_BJYY.append(k+" "+i)
				elif k in io_gene_N:
					if not re.search("融合", i):
						io_n_list_BJYY.append(k+" "+i)
					else:
						io_n_list_BJYY.append(i)
	if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF19" in io_result.keys():
		io_n_list_BJYY.append("CCND1/FGF3/FGF19扩增")
	io_p_result = ", ".join(io_p_list_BJYY)
	io_n_result = ", ".join(io_n_list_BJYY)
	num = len(io_p_list_BJYY+io_n_list_BJYY)
	
	return io_result, io_p_result, io_n_result, num
# 合并结束-2023.11.20

# 2026.03.17-新增一个格式用于北大人民MP
def io_detect_for_bdrm(var_data):
	# io_result用于填充IO表
	io_result = {}
	# io_p_result用于填充检测结果小结
	io_p_result = ""
	# io_n_result用于填充检测结果小结
	io_n_result = ""
	# 汇总体细胞I/II/肿瘤发生发展相关变异+胚系致病/疑似致病性变异
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","ERCC1","FANCA","MRE11","PALB2","RAD50","XRCC1","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53","KRAS","CD274","ARID1A","LRP1B","SETD2","PRKDC","TERT","KMT2D","FAT1","CDK12","SERPINB3","SERPINB4"]
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1","IRF1","JAK1","JAK2","APC","CTNNB1","B2M","PTEN","CCND1","FGF3","FGF19"]
	level_12_var = [var for var in var_data if (var["var_origin"] != "germline" and var["clinic_num_s"] in [5, 4]) or (var["var_origin"] == "germline" and var["clinic_num_g"] in [5, 4])]
	# 2025.02.21- 更新展示变异等级
	# 1. 有配对：体细胞I/II/肿瘤发生发展 + 胚系4/5 + 胚系3类但有用药
	# 2. 无配对：体细胞I/II/III类 + 预测胚系4/5（归为肿瘤发生发展） + 预测胚系3类但有用药
	# ==> 即新增（确认/预测）胚系3类但有用药的变异即可
	germline_level3_regimen = [var for var in var_data if var["var_origin"] == "germline" and var["clinic_num_g"] == 3 and judgeRegimen(var)]
	if germline_level3_regimen:
		level_12_var.extend(germline_level3_regimen)
	# 2025.02.21-新增完成

	for var in level_12_var:
		# 仅展示扩增的基因
		# 扩增共突变
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in ["CD274", "MDM2", "MDM4", "CCND1", "FGF3", "FGF19"]:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append("扩增")
		# 仅展示融合的基因
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			io_result["ALK"].append(var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"融合")
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["hgvs_p"] != "p.?":
				io_result[var["gene_symbol"]].append(var["hgvs_p_BJYY"])
			else:
				io_result[var["gene_symbol"]].append(var["hgvs_c"])
	
	# summary展示
	io_p_list_BJYY = []
	io_n_list_BJYY = []
	for k, v in io_result.items():
		if k not in ["CCND1","FGF3","FGF19"]:
			for i in v:
				if k in io_gene_P:
					io_p_list_BJYY.append(k+" : "+i)
				elif k in io_gene_N:
					if not re.search("融合", i):
						io_n_list_BJYY.append(k+" : "+i)
					else:
						io_n_list_BJYY.append(i)
	if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF19" in io_result.keys():
		io_n_list_BJYY.append("CCND1/FGF3/FGF19扩增")
	io_p_result = ", ".join(io_p_list_BJYY)
	io_n_result = ", ".join(io_n_list_BJYY)
	num = len(io_p_list_BJYY+io_n_list_BJYY)
	
	return io_result, io_p_result, io_n_result, num
# 2026.03.17-新增完成

# 北京医院-合并梦晨代码-2023.11.20
# 2025.07.21-增加HD-刘炜芬
def io_detect_for_BJYY_hd(var_data):
	# io_result用于填充IO表
	io_result = {}
	# io_p_result用于填充检测结果小结
	io_p_result = ""
	# io_n_result用于填充检测结果小结
	io_n_result = ""
	# 汇总体细胞I/II/肿瘤发生发展相关变异+胚系致病/疑似致病性变异
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","ERCC1","FANCA","MRE11","PALB2","RAD50","XRCC1","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53","KRAS","CD274","ARID1A","LRP1B","SETD2","PRKDC","TERT","KMT2D","FAT1","CDK12","SERPINB3","SERPINB4"]
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1","IRF1","JAK1","JAK2","APC","CTNNB1","B2M","PTEN","CCND1","FGF3","FGF19"]
	# hd 是额外需要展示的，原有的输出类型不变
	hd_gene_list = ["ATM", "BRCA1", "BRCA2", "BRIP1", "CDK12", "CHEK1", "CHEK2", "FANCA", \
				 	"PALB2", "SETD2", "TP53", "CDKN2A", "CDKN2B", "PTEN", "STK11"]
	level_12_var = [var for var in var_data if (var["var_origin"] != "germline" and var["clinic_num_s"] in [5, 4]) or (var["var_origin"] == "germline" and var["clinic_num_g"] in [5, 4])]
	# 2025.02.21- 更新展示变异等级
	# 1. 有配对：体细胞I/II/肿瘤发生发展 + 胚系4/5 + 胚系3类但有用药
	# 2. 无配对：体细胞I/II/III类 + 预测胚系4/5（归为肿瘤发生发展） + 预测胚系3类但有用药
	# ==> 即新增（确认/预测）胚系3类但有用药的变异即可
	germline_level3_regimen = [var for var in var_data if var["var_origin"] == "germline" and var["clinic_num_g"] == 3 and judgeRegimen(var)]
	if germline_level3_regimen:
		level_12_var.extend(germline_level3_regimen)
	# 2025.02.21-新增完成

	num = 0
	for var in level_12_var:
		# 仅展示扩增的基因
		# 扩增共突变
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in ["CD274", "MDM2", "MDM4", "CCND1", "FGF3", "FGF19"]:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append("扩增")
			if var["gene_symbol"] in ["CD274", "MDM2", "MDM4"]:
				num += 1
		# 仅展示融合的基因
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			io_result["ALK"].append(var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"融合")
			num += 1
		# HD基因经确认同时展示snvindel和hd
		elif (var["bio_category"] == "Snvindel" or var["bio_category"] == "PHd") and var["gene_symbol"] in hd_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["bio_category"] == "Snvindel":
				if var["hgvs_p"] != "p.?":
					io_result[var["gene_symbol"]].append(var["hgvs_p_BJYY"])
				else:
					io_result[var["gene_symbol"]].append(var["hgvs_c"])
			elif var["bio_category"] == "PHd":
				if var["type"] == "HomoDel":
					if "纯合缺失" not in io_result[var["gene_symbol"]]:
						io_result[var["gene_symbol"]].append("纯合缺失")
				elif var["type"] == "HeteDel":
					if "杂合缺失" not in io_result[var["gene_symbol"]]:
						io_result[var["gene_symbol"]].append("杂合缺失")
				else:
					io_result[var["gene_symbol"]].append("未知变异类型！")
			num += 1
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["hgvs_p"] != "p.?":
				io_result[var["gene_symbol"]].append(var["hgvs_p_BJYY"])
			else:
				io_result[var["gene_symbol"]].append(var["hgvs_c"])
			num += 1
	
	# summary展示
	io_p_list_BJYY = []
	io_n_list_BJYY = []
	for k, v in io_result.items():
		if k not in ["CCND1","FGF3","FGF19"]:
			for i in v:
				if k in io_gene_P:
					io_p_list_BJYY.append(k+" "+i)
				elif k in io_gene_N:
					if not re.search("融合", i):
						io_n_list_BJYY.append(k+" "+i)
					else:
						io_n_list_BJYY.append(i)
	if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF19" in io_result.keys():
		io_n_list_BJYY.append("CCND1/FGF3/FGF19扩增")
		num += 1
	io_p_result = ", ".join(io_p_list_BJYY)
	io_n_result = ", ".join(io_n_list_BJYY)

	# 2025.07.21-HD同一个基因可能返回多个变异，去重后统计不准
	#num = len(io_p_list_BJYY+io_n_list_BJYY)
	
	return io_result, io_p_result, io_n_result, num
# 合并结束-2023.11.20

# 孙逸仙116-三字母-2024.08.19
def io_detect_for_116_syx(var_data):
	# 返回结果中的io_result用于填充IO表，", ".join(io_p_list), ", ".join(io_n_list)用于填充检测结果小结
	io_result = {}
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","FANCA","MRE11","PALB2","MLH1","MSH2","MSH6",\
				 "PMS2","POLE","POLD1","TP53","KRAS","CD274","ARID1A","TERT","CDK12"]
	#io_gene_N = ["EGFR","ALK","CDKN2A","CDKN2B","STK11","JAK1","JAK2","APC","CTNNB1","PTEN"]
	# 新增CCND1/FGF3/FGF19共扩增-2023.07.12
	io_gene_N = ["EGFR","ALK","CDKN2A","CDKN2B","STK11","JAK1","JAK2","APC","CTNNB1","PTEN", "CCND1", "FGF3", "FGF19"]
	level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]
	# 2025.02.21- 更新展示变异等级
	# 1. 有配对：体细胞I/II/肿瘤发生发展 + 胚系4/5 + 胚系3类但有用药
	# 2. 无配对：体细胞I/II/肿瘤发生发展 + 预测胚系4/5（归为I/II/肿瘤发生发展） + 预测胚系3类但有用药
	# ==> 即新增（确认/预测）胚系3类但有用药的变异即可
	germline_level3_regimen = [var for var in var_data if var["var_origin"] == "germline" and var["clinic_num_g"] == 3 and judgeRegimen(var)]
	if germline_level3_regimen:
		level_12_var.extend(germline_level3_regimen)
	# 2025.02.21-新增完成
	for var in level_12_var:
		# 仅展示扩增的基因
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in ["CD274", "CCND1", "FGF3", "FGF19"]:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append("扩增")
		# 仅展示融合的基因
		if var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			io_result["ALK"].append(var["five_prime_gene"]+"-"+var["three_prime_gene"]+"融合")
			# 加一个孙逸仙的 Gene1-ALK(G8:A12)基因重排 -2023.08.02
			if "ALK_syx" not in io_result.keys():
				io_result.setdefault("ALK_syx", [])
			io_result["ALK_syx"].append("{0}-{1}({2}{3}:{4}{5})基因重排".format(var["five_prime_gene"], var["three_prime_gene"], var["five_prime_gene"][0], \
								var["five_prime_cds"].replace("exon", ""), var["three_prime_gene"][0], var["three_prime_cds"].replace("exon", "")))
			# 新增结束-2023.08.02
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["hgvs_p_abbr"] and var["hgvs_p_abbr"] != "p.?":
				io_result[var["gene_symbol"]].append(var["hgvs_p_abbr"].replace("(", "").replace(")", ""))
			else:
				io_result[var["gene_symbol"]].append(var["hgvs_c"])
	# summary展示
	both_cnv_list = ["CCND1", "FGF3", "FGF19"]
	io_p_list = ["{0} {1}".format(k, i) for k,v in io_result.items() for i in v if k in io_gene_P]
	io_n_list = ["{0} {1}".format(k, i) if not re.search("融合", i) else i for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_N]
	# 处理CNV共突变
	if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF19" in io_result.keys():
		io_n_list.append("CCND1/FGF3/FGF19扩增")
	
	return io_result, ", ".join(io_p_list), ", ".join(io_n_list)

def io_detect_for_cp200(var_data):
	# 返回结果中的io_result用于填充IO表，", ".join(io_p_list), ", ".join(io_n_list)用于填充检测结果小结
	io_result = {}
	# 汇总体细胞或来源不明I/II/肿瘤发生发展相关变异
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","FANCA","MRE11",\
				 "PALB2","RAD50","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53",\
				 "KRAS","CD274","ARID1A","SETD2","TERT","KMT2D","FAT1","CDK12"]
	
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1",\
				 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","FGF19"]
	cnv_gene_list = ["CD274", "MDM2", "MDM4", "FGF19"]
	
	level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]

	for var in level_12_var:
		# 仅展示扩增的基因
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in cnv_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append("扩增")
		# 仅展示融合的基因
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			if var["five_prime_gene"]+"-"+var["three_prime_gene"]+"融合" not in io_result["ALK"]:
				io_result["ALK"].append(var["five_prime_gene"]+"-"+var["three_prime_gene"]+"融合")
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["hgvs_p"] != "p.?":
				io_result[var["gene_symbol"]].append(var["hgvs_p"])
			else:
				io_result[var["gene_symbol"]].append(var["hgvs_c"])
	# summary展示
	io_p_list = ["{0} {1}".format(k, i) for k,v in io_result.items() for i in v if k in io_gene_P]
	io_n_list = ["{0} {1}".format(k, i) if not re.search("融合", i) else i for k,v in io_result.items() for i in v if k in io_gene_N]

	return ", ".join(io_p_list), ", ".join(io_n_list)

def io_detect_for_cp200_zsly(var_data):
	# 免疫负相关删除了CDKN2B
	# 返回结果中的io_result用于填充IO表，", ".join(io_p_list), ", ".join(io_n_list)用于填充检测结果小结
	io_result = {}
	# 汇总体细胞或来源不明I/II/肿瘤发生发展相关变异
	#io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","FANCA","MRE11",\
	#			 "PALB2","RAD50","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53",\
	#			 "KRAS","CD274","ARID1A","SETD2","TERT","KMT2D","FAT1","CDK12"]
	#2026.05.07 增加DNMT3A 孟智悦
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","FANCA","MRE11",\
				 "PALB2","RAD50","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53",\
				 "KRAS","CD274","ARID1A","SETD2","TERT","KMT2D","FAT1","CDK12","DNMT3A"]

	#io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","DNMT3A","STK11","IFNGR1",\
	#			 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","FGF19"]
	# 2026.05.07 删除基因DNMT3A 孟智悦
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","STK11","IFNGR1",\
				 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","FGF19"]

	cnv_gene_list = ["CD274", "MDM2", "MDM4", "FGF19"]
	
	level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]

	for var in level_12_var:
		# 仅展示扩增的基因
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in cnv_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append("扩增")
		# 仅展示融合的基因
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			if var["five_prime_gene"]+"-"+var["three_prime_gene"]+"融合" not in io_result["ALK"]:
				io_result["ALK"].append(var["five_prime_gene"]+"-"+var["three_prime_gene"]+"融合")
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["hgvs_p"] != "p.?":
				io_result[var["gene_symbol"]].append(var["hgvs_p"])
			else:
				io_result[var["gene_symbol"]].append(var["hgvs_c"])
	# summary展示
	io_p_list = ["{0} {1}".format(k, i) for k,v in io_result.items() for i in v if k in io_gene_P]
	io_n_list = ["{0} {1}".format(k, i) if not re.search("融合", i) else i for k,v in io_result.items() for i in v if k in io_gene_N]

	return ", ".join(io_p_list), ", ".join(io_n_list)

def io_detect_for_boncopro(var_data):
	# 返回结果中的io_result用于填充IO表，", ".join(io_p_list), ", ".join(io_n_list)用于填充检测结果小结
	io_result = {}
	# 汇总体细胞I/II/肿瘤发生发展相关变异+胚系致病/疑似致病性变异
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","ERCC1","FANCA","MRE11",\
				 "PALB2","RAD50","XRCC1","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53",\
				 "KRAS","CD274","ARID1A","SETD2","TERT","CDK12"]
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","STK11",\
				 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN"]
	cnv_gene_list = ["CD274", "MDM2", "MDM4"]
	
	level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]
	# 2025.02.21- 更新展示变异等级
	# 1. 有配对：体细胞I/II/肿瘤发生发展 + 胚系4/5 + 胚系3类但有用药
	# 2. 无配对：体细胞I/II/肿瘤发生发展 + 预测胚系4/5（归为I/II/肿瘤发生发展） + 预测胚系3类但有用药
	# ==> 即新增（确认/预测）胚系3类但有用药的变异即可
	germline_level3_regimen = [var for var in var_data if var["var_origin"] == "germline" and var["clinic_num_g"] == 3 and judgeRegimen(var)]
	if germline_level3_regimen:
		level_12_var.extend(germline_level3_regimen)
	# 2025.02.21-新增完成

	for var in level_12_var:
		# 仅展示扩增的基因
		# 2026.01.12-CNV 需要排除缺失的情况
		#if var["bio_category"] == "Cnv" and var["gene_symbol"] in cnv_gene_list:
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in cnv_gene_list and "cnv_type" in var.keys() and var["cnv_type"] and var["cnv_type"] not in ["Loss", "loss"]:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append("扩增")
		# 仅展示融合的基因
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			if var["five_prime_gene"]+"-"+var["three_prime_gene"]+"融合" not in io_result["ALK"]:
				io_result["ALK"].append(var["five_prime_gene"]+"-"+var["three_prime_gene"]+"融合")
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["hgvs_p"] != "p.?":
				io_result[var["gene_symbol"]].append(var["hgvs_p"])
			else:
				io_result[var["gene_symbol"]].append(var["hgvs_c"])
	# summary展示
	io_p_list = ["{0} {1}".format(k, i) for k,v in io_result.items() for i in v if k in io_gene_P]
	io_n_list = ["{0} {1}".format(k, i) if not re.search("融合", i) else i for k,v in io_result.items() for i in v if k in io_gene_N]

	# 新增返回io_result-2026.01.12
	#return ", ".join(io_p_list), ", ".join(io_n_list)
	return io_result, ", ".join(io_p_list), ", ".join(io_n_list)

# 中山六院io结果需要频率/拷贝数-2025.04.02
def ZSLY_io_detect(var_data):
	io_result = {}
	# 汇总体细胞I/II/肿瘤发生发展相关变异
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","FANCA","MRE11",\
				 "PALB2","RAD50","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53",\
				 "KRAS","CD274","ARID1A","SETD2","TERT","KMT2D","FAT1","CDK12"]
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1",\
				 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","FGF19"]
	cnv_gene_list = ["CD274", "MDM2", "MDM4", "FGF19"]

	level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]
	
	for var in level_12_var:
		# 仅展示扩增的基因
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in cnv_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append(var)
		# 仅展示融合的基因
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			io_result["ALK"].append(var)
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append(var)

	return io_result

# 2025.06.17-新增，包含HD结果
def io_detect_hd(var_data):
	# 返回结果中的io_result用于填充IO表，", ".join(io_p_list), ", ".join(io_n_list)用于填充检测结果小结
	io_result = {}
	# 汇总体细胞I/II/肿瘤发生发展相关变异+胚系致病/疑似致病性变异
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","ERCC1","FANCA","MRE11",\
				 "PALB2","RAD50","XRCC1","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53",\
				 "KRAS","CD274","ARID1A","LRP1B","SETD2","PRKDC","TERT","KMT2D","FAT1","CDK12",\
				 "SERPINB3","SERPINB4"]
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1","IRF1",\
				 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","CCND1","FGF3","FGF19"]
	cnv_gene_list = ["CD274", "MDM2", "MDM4", "CCND1", "FGF3", "FGF19"]
	# hd 是额外需要展示的，原有的输出类型不变
	hd_gene_list = ["ATM", "BRCA1", "BRCA2", "BRIP1", "CDK12", "CHEK1", "CHEK2", "FANCA", \
				 	"PALB2", "SETD2", "TP53", "CDKN2A", "CDKN2B", "PTEN", "STK11"]

	level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]
	# 2025.02.21- 更新展示变异等级
	# 1. 有配对：体细胞I/II/肿瘤发生发展 + 胚系4/5 + 胚系3类但有用药
	# 2. 无配对：体细胞I/II/肿瘤发生发展 + 预测胚系4/5（归为I/II/肿瘤发生发展） + 预测胚系3类但有用药
	# ==> 即新增（确认/预测）胚系3类但有用药的变异即可
	germline_level3_regimen = [var for var in var_data if var["var_origin"] == "germline" and var["clinic_num_g"] == 3 and judgeRegimen(var)]
	if germline_level3_regimen:
		level_12_var.extend(germline_level3_regimen)
	# 2025.02.21-新增完成

	for var in level_12_var:
		# 仅展示扩增的基因
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in cnv_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append("扩增")
		# 仅展示融合的基因
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			# 融合可能会出现exon相同断点不同的变异，报告中会重复，这边加个去重-2023.04.13
			if var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"融合" not in io_result["ALK"]:
				io_result["ALK"].append(var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"融合")
		# HD基因经确认同时展示snvindel和hd
		elif (var["bio_category"] == "Snvindel" or var["bio_category"] == "PHd") and var["gene_symbol"] in hd_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["bio_category"] == "Snvindel":
				if var["hgvs_p"] != "p.?":
					io_result[var["gene_symbol"]].append(var["hgvs_p"])
				else:
					io_result[var["gene_symbol"]].append(var["hgvs_c"])
			elif var["bio_category"] == "PHd":
				if var["type"] == "HomoDel":
					if "纯合缺失" not in io_result[var["gene_symbol"]]:
						io_result[var["gene_symbol"]].append("纯合缺失")
				elif var["type"] == "HeteDel":
					if "杂合缺失" not in io_result[var["gene_symbol"]]:
						io_result[var["gene_symbol"]].append("杂合缺失")
				else:
					io_result[var["gene_symbol"]].append("未知变异类型！")
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["hgvs_p"] != "p.?":
				io_result[var["gene_symbol"]].append(var["hgvs_p"])
			else:
				io_result[var["gene_symbol"]].append(var["hgvs_c"])
	# summary展示
	both_cnv_list = ["CCND1","FGF3","FGF19"]
	io_p_list = ["{0} {1}".format(k, i) for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_P]
	io_n_list = ["{0} {1}".format(k, i) if not re.search("融合", i) else i for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_N]

	# 处理CNV共突变
	if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF19" in io_result.keys():
		io_n_list.append("CCND1/FGF3/FGF19共扩增")

	return io_result, ", ".join(io_p_list), ", ".join(io_n_list)

# 2025.07.17-增加HD，其他不变
# 2026.03.19-CNV/SV/HD增加转录本
def io_detect_for_ZJZL_hd(var_data, config):
	fdzs_dict, fjzl_database, Data = getconfigxlsx(config)
	# io_result用于填充IO表
	io_result = {}

	# 删除FGF4扩增-2025.02.21
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","ERCC1","FANCA","MRE11",\
				 "PALB2","RAD50","XRCC1","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53",\
				 "KRAS","CD274","ARID1A","LRP1B","SETD2","PRKDC","TERT","KMT2D","FAT1","CDK12",\
				 "SERPINB3","SERPINB4"]
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1","IRF1",\
				 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","CCND1","FGF3","FGF19"]
	cnv_gene_list = ["CD274", "MDM2", "MDM4", "CCND1", "FGF3", "FGF19"]
	# hd 是额外需要展示的，原有的输出类型不变
	hd_gene_list = ["ATM", "BRCA1", "BRCA2", "BRIP1", "CDK12", "CHEK1", "CHEK2", "FANCA", \
				 	"PALB2", "SETD2", "TP53", "CDKN2A", "CDKN2B", "PTEN", "STK11"]
	# 展示I/II类和肿瘤发生发展相关变异
	level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]
	# 2025.02.21- 更新展示变异等级
	# 1. 有配对：体细胞I/II/肿瘤发生发展 + 胚系4/5 + 胚系3类但有用药
	# 2. 无配对：体细胞I/II/肿瘤发生发展 + 预测胚系4/5（归为I/II/肿瘤发生发展） + 预测胚系3类但有用药
	# ==> 即新增（确认/预测）胚系3类但有用药的变异即可
	germline_level3_regimen = [var for var in var_data if var["var_origin"] == "germline" and var["clinic_num_g"] == 3 and judgeRegimen(var)]
	if germline_level3_regimen:
		level_12_var.extend(germline_level3_regimen)
	# 2025.02.21-新增完成
	for var in level_12_var:
		# 仅展示扩增的基因
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in cnv_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append(
				{
					"var_type" : "cnv", 
					"var_info" : "扩增",
					"transcript_primary" : var["transcript_primary"] if "transcript_primary" in var.keys() and var["transcript_primary"] else ""
					}
				)
		# 仅展示融合的基因
		# 可能存在两个融合基因都在检测范围里的情况，系统json返回的gene_symbol为gene1,gene2，需额外做识别
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			io_result["ALK"].append(
				{
					"var_type" : "sv", 
					"three_prime_gene" : var["three_prime_gene"], 
					"three_prime_cds" : var["three_prime_cds"], 
					"five_prime_gene" : var["five_prime_gene"], 
					"five_prime_cds" : var["five_prime_cds"],
					"five_prime_transcript" : var["five_prime_transcript"],
					"three_prime_transcript" : var["three_prime_transcript"]
					}
				)
		# HD基因经确认同时展示snvindel和hd
		elif (var["bio_category"] == "Snvindel" or var["bio_category"] == "PHd") and var["gene_symbol"] in hd_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["bio_category"] == "Snvindel":
				io_result[var["gene_symbol"]].append(
				{
					"var_type" : "snvindel", 
					"hgvs_p" : var["hgvs_p"], 
					"hgvs_c" : var["hgvs_c"], 
					"var_origin" : var["var_origin"], 
					"gene_region" : var["gene_region"], 
					"transcript_primary" : var["transcript_primary"]
					}
				)
			elif var["bio_category"] == "PHd":
				io_result[var["gene_symbol"]].append(
					{
						"var_type" : "PHd",
						"type" : var["type"],
						"region" : var["region"],
						"transcript_primary" : var["transcript_primary"]
					}
				)
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append(
				{
					"var_type" : "snvindel", 
					"hgvs_p" : var["hgvs_p"], 
					"hgvs_c" : var["hgvs_c"], 
					"var_origin" : var["var_origin"], 
					"gene_region" : var["gene_region"], 
					"transcript_primary" : var["transcript_primary"]
					}
				)
	# summary展示
	io_inter = []
	for k, v in io_result.items():
		if k not in ["CCND1","FGF3","FGF19"]:
			if k in Data.keys():
				io_inter.append({
					"gene_symbol" : k,
					"var_info" : v,
					"inter" : Data[k]
				})
	if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF19" in io_result.keys():
		io_inter.append({
			"gene_symbol" : "CCND1/FGF3/FGF19",
			"var_info" : io_result["CCND1"],
			"inter" : Data["CCND1"] 
		})
	
	io_p_list = [i for i in io_inter if i["gene_symbol"] in io_gene_P]
	io_n_list = [i for i in io_inter if i["gene_symbol"] in io_gene_N or i["gene_symbol"] == "CCND1/FGF3/FGF19"]

	return io_p_list, io_n_list

# 2026.02.03-适用武汉协和MP，氨基酸要展示三字母缩写
def abbr_io_detect_hd(var_data):
	# 返回结果中的io_result用于填充IO表，", ".join(io_p_list), ", ".join(io_n_list)用于填充检测结果小结
	io_result = {}
	# 汇总体细胞I/II/肿瘤发生发展相关变异+胚系致病/疑似致病性变异
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","ERCC1","FANCA","MRE11",\
				 "PALB2","RAD50","XRCC1","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53",\
				 "KRAS","CD274","ARID1A","LRP1B","SETD2","PRKDC","TERT","KMT2D","FAT1","CDK12",\
				 "SERPINB3","SERPINB4"]
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1","IRF1",\
				 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","CCND1","FGF3","FGF19"]
	cnv_gene_list = ["CD274", "MDM2", "MDM4", "CCND1", "FGF3", "FGF19"]
	# hd 是额外需要展示的，原有的输出类型不变
	hd_gene_list = ["ATM", "BRCA1", "BRCA2", "BRIP1", "CDK12", "CHEK1", "CHEK2", "FANCA", \
				 	"PALB2", "SETD2", "TP53", "CDKN2A", "CDKN2B", "PTEN", "STK11"]

	level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]
	# 2025.02.21- 更新展示变异等级
	# 1. 有配对：体细胞I/II/肿瘤发生发展 + 胚系4/5 + 胚系3类但有用药
	# 2. 无配对：体细胞I/II/肿瘤发生发展 + 预测胚系4/5（归为I/II/肿瘤发生发展） + 预测胚系3类但有用药
	# ==> 即新增（确认/预测）胚系3类但有用药的变异即可
	germline_level3_regimen = [var for var in var_data if var["var_origin"] == "germline" and var["clinic_num_g"] == 3 and judgeRegimen(var)]
	if germline_level3_regimen:
		level_12_var.extend(germline_level3_regimen)
	# 2025.02.21-新增完成

	for var in level_12_var:
		# 仅展示扩增的基因
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in cnv_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append("扩增")
		# 仅展示融合的基因
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			# 融合可能会出现exon相同断点不同的变异，报告中会重复，这边加个去重-2023.04.13
			if var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"融合" not in io_result["ALK"]:
				io_result["ALK"].append(var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"融合")
		# HD基因经确认同时展示snvindel和hd
		elif (var["bio_category"] == "Snvindel" or var["bio_category"] == "PHd") and var["gene_symbol"] in hd_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["bio_category"] == "Snvindel":
				if var["hgvs_p"] != "p.?":
					io_result[var["gene_symbol"]].append(var["hgvs_p_abbr"])
				else:
					io_result[var["gene_symbol"]].append(var["hgvs_c"])
			elif var["bio_category"] == "PHd":
				if var["type"] == "HomoDel":
					if "纯合缺失" not in io_result[var["gene_symbol"]]:
						io_result[var["gene_symbol"]].append("纯合缺失")
				elif var["type"] == "HeteDel":
					if "杂合缺失" not in io_result[var["gene_symbol"]]:
						io_result[var["gene_symbol"]].append("杂合缺失")
				else:
					io_result[var["gene_symbol"]].append("未知变异类型！")
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["hgvs_p"] != "p.?":
				io_result[var["gene_symbol"]].append(var["hgvs_p_abbr"])
			else:
				io_result[var["gene_symbol"]].append(var["hgvs_c"])
	# summary展示
	both_cnv_list = ["CCND1","FGF3","FGF19"]
	io_p_list = ["{0} {1}".format(k, i) for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_P]
	io_n_list = ["{0} {1}".format(k, i) if not re.search("融合", i) else i for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_N]

	# 处理CNV共突变
	if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF19" in io_result.keys():
		io_n_list.append("CCND1/FGF3/FGF19共扩增")

	return io_result, ", ".join(io_p_list), ", ".join(io_n_list)


# 2026.02.25-新增，包含HD结果，适用北大三MP，胚系3类有用药的位点不展示
def io_detect_hd_bds(var_data):
	# 返回结果中的io_result用于填充IO表，", ".join(io_p_list), ", ".join(io_n_list)用于填充检测结果小结
	io_result = {}
	# 汇总体细胞I/II/肿瘤发生发展相关变异+胚系致病/疑似致病性变异
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","ERCC1","FANCA","MRE11",\
				 "PALB2","RAD50","XRCC1","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53",\
				 "KRAS","CD274","ARID1A","LRP1B","SETD2","PRKDC","TERT","KMT2D","FAT1","CDK12",\
				 "SERPINB3","SERPINB4"]
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1","IRF1",\
				 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","CCND1","FGF3","FGF19"]
	cnv_gene_list = ["CD274", "MDM2", "MDM4", "CCND1", "FGF3", "FGF19"]
	# hd 是额外需要展示的，原有的输出类型不变
	hd_gene_list = ["ATM", "BRCA1", "BRCA2", "BRIP1", "CDK12", "CHEK1", "CHEK2", "FANCA", \
				 	"PALB2", "SETD2", "TP53", "CDKN2A", "CDKN2B", "PTEN", "STK11"]

	level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]
	# 2025.02.21- 更新展示变异等级
	# 1. 有配对：体细胞I/II/肿瘤发生发展 + 胚系4/5 + 胚系3类但有用药
	# 2. 无配对：体细胞I/II/肿瘤发生发展 + 预测胚系4/5（归为I/II/肿瘤发生发展） + 预测胚系3类但有用药
	# ==> 即新增（确认/预测）胚系3类但有用药的变异即可
	#germline_level3_regimen = [var for var in var_data if var["var_origin"] == "germline" and var["clinic_num_g"] == 3 and judgeRegimen(var)]
	#if germline_level3_regimen:
	#	level_12_var.extend(germline_level3_regimen)
	# 2025.02.21-新增完成

	for var in level_12_var:
		# 仅展示扩增的基因
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in cnv_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append("扩增")
		# 仅展示融合的基因
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			# 融合可能会出现exon相同断点不同的变异，报告中会重复，这边加个去重-2023.04.13
			if var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"融合" not in io_result["ALK"]:
				io_result["ALK"].append(var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"融合")
		# HD基因经确认同时展示snvindel和hd
		elif (var["bio_category"] == "Snvindel" or var["bio_category"] == "PHd") and var["gene_symbol"] in hd_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["bio_category"] == "Snvindel":
				if var["hgvs_p"] != "p.?":
					io_result[var["gene_symbol"]].append(var["hgvs_p"])
				else:
					io_result[var["gene_symbol"]].append(var["hgvs_c"])
			elif var["bio_category"] == "PHd":
				if var["type"] == "HomoDel":
					if "纯合缺失" not in io_result[var["gene_symbol"]]:
						io_result[var["gene_symbol"]].append("纯合缺失")
				elif var["type"] == "HeteDel":
					if "杂合缺失" not in io_result[var["gene_symbol"]]:
						io_result[var["gene_symbol"]].append("杂合缺失")
				else:
					io_result[var["gene_symbol"]].append("未知变异类型！")
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["hgvs_p"] != "p.?":
				io_result[var["gene_symbol"]].append(var["hgvs_p"])
			else:
				io_result[var["gene_symbol"]].append(var["hgvs_c"])
	# summary展示
	both_cnv_list = ["CCND1","FGF3","FGF19"]
	io_p_list = ["{0} {1}".format(k, i) for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_P]
	io_n_list = ["{0} {1}".format(k, i) if not re.search("融合", i) else i for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_N]

	# 处理CNV共突变
	if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF19" in io_result.keys():
		io_n_list.append("CCND1/FGF3/FGF19共扩增")

	return io_result, ", ".join(io_p_list), ", ".join(io_n_list)

# 2026.05.08 新增，北大人民，包含HD结果
def io_detect_for_bdrm_hd(var_data):
	# io_result用于填充IO表
	io_result = {}
	# io_p_result用于填充检测结果小结
	io_p_result = ""
	# io_n_result用于填充检测结果小结
	io_n_result = ""
	# 汇总体细胞I/II/肿瘤发生发展相关变异+胚系致病/疑似致病性变异
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","ERCC1","FANCA","MRE11","PALB2","RAD50","XRCC1","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53","KRAS","CD274","ARID1A","LRP1B","SETD2","PRKDC","TERT","KMT2D","FAT1","CDK12","SERPINB3","SERPINB4"]
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1","IRF1","JAK1","JAK2","APC","CTNNB1","B2M","PTEN","CCND1","FGF3","FGF19"]
	# hd 是额外需要展示的，原有的输出类型不变
	hd_gene_list = ["ATM", "BRCA1", "BRCA2", "BRIP1", "CDK12", "CHEK1", "CHEK2", "FANCA", \
				 	"PALB2", "SETD2", "TP53", "CDKN2A", "CDKN2B", "PTEN", "STK11"]
	level_12_var = [var for var in var_data if (var["var_origin"] != "germline" and var["clinic_num_s"] in [5, 4]) or (var["var_origin"] == "germline" and var["clinic_num_g"] in [5, 4])]
	# 2025.02.21- 更新展示变异等级
	# 1. 有配对：体细胞I/II/肿瘤发生发展 + 胚系4/5 + 胚系3类但有用药
	# 2. 无配对：体细胞I/II/III类 + 预测胚系4/5（归为肿瘤发生发展） + 预测胚系3类但有用药
	# ==> 即新增（确认/预测）胚系3类但有用药的变异即可
	germline_level3_regimen = [var for var in var_data if var["var_origin"] == "germline" and var["clinic_num_g"] == 3 and judgeRegimen(var)]
	if germline_level3_regimen:
		level_12_var.extend(germline_level3_regimen)
	# 2025.02.21-新增完成

	num = 0
	for var in level_12_var:
		# 仅展示扩增的基因
		# 扩增共突变
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in ["CD274", "MDM2", "MDM4", "CCND1", "FGF3", "FGF19"]:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append("扩增")
			if var["gene_symbol"] in ["CD274", "MDM2", "MDM4"]:
				num += 1
		# 仅展示融合的基因
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			io_result["ALK"].append(var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"融合")
			num += 1
		# HD基因经确认同时展示snvindel和hd
		elif (var["bio_category"] == "Snvindel" or var["bio_category"] == "PHd") and var["gene_symbol"] in hd_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["bio_category"] == "Snvindel":
				if var["hgvs_p"] != "p.?":
					io_result[var["gene_symbol"]].append(var["hgvs_p_BJYY"])
				else:
					io_result[var["gene_symbol"]].append(var["hgvs_c"])
			elif var["bio_category"] == "PHd":
				if var["type"] == "HomoDel":
					if "纯合缺失" not in io_result[var["gene_symbol"]]:
						io_result[var["gene_symbol"]].append("纯合缺失")
				elif var["type"] == "HeteDel":
					if "杂合缺失" not in io_result[var["gene_symbol"]]:
						io_result[var["gene_symbol"]].append("杂合缺失")
				else:
					io_result[var["gene_symbol"]].append("未知变异类型！")
			num += 1
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["hgvs_p"] != "p.?":
				io_result[var["gene_symbol"]].append(var["hgvs_p_BJYY"])
			else:
				io_result[var["gene_symbol"]].append(var["hgvs_c"])
			num += 1
	
	# summary展示
	io_p_list_BJYY = []
	io_n_list_BJYY = []
	for k, v in io_result.items():
		if k not in ["CCND1","FGF3","FGF19"]:
			for i in v:
				if k in io_gene_P:
					io_p_list_BJYY.append(k+" : "+i)
				elif k in io_gene_N:
					if not re.search("融合", i):
						io_n_list_BJYY.append(k+" : "+i)
					else:
						io_n_list_BJYY.append(i)
	if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF19" in io_result.keys():
		io_n_list_BJYY.append("CCND1/FGF3/FGF19扩增")
		num += 1
	io_p_result = ", ".join(io_p_list_BJYY)
	io_n_result = ", ".join(io_n_list_BJYY)
	#num = len(io_p_list_BJYY+io_n_list_BJYY)
	
	return io_result, io_p_result, io_n_result, num

# 复旦中山mp，cnv描述修改，嵇梦晨，2026.05.14
def io_detect_fdzs_mp(var_data):
	# 返回结果中的io_result用于填充IO表，", ".join(io_p_list), ", ".join(io_n_list)用于填充检测结果小结
	io_result = {}
	# 汇总体细胞I/II/肿瘤发生发展相关变异+胚系致病/疑似致病性变异
	io_gene_P = ["ATM","ATR","BRCA1","BRCA2","BRIP1","CHEK1","CHEK2","ERCC1","FANCA","MRE11",\
				 "PALB2","RAD50","XRCC1","MLH1","MSH2","MSH6","PMS2","POLE","POLD1","TP53",\
				 "KRAS","CD274","ARID1A","LRP1B","SETD2","PRKDC","TERT","KMT2D","FAT1","CDK12",\
				 "SERPINB3","SERPINB4"]
	io_gene_N = ["EGFR","ALK","MDM2","MDM4","CDKN2A","CDKN2B","DNMT3A","STK11","IFNGR1","IRF1",\
				 "JAK1","JAK2","APC","CTNNB1","B2M","PTEN","CCND1","FGF3","FGF19"]
	cnv_gene_list = ["CD274", "MDM2", "MDM4", "CCND1", "FGF3", "FGF19"]
	# hd 是额外需要展示的，原有的输出类型不变
	hd_gene_list = ["ATM", "BRCA1", "BRCA2", "BRIP1", "CDK12", "CHEK1", "CHEK2", "FANCA", \
				 	"PALB2", "SETD2", "TP53", "CDKN2A", "CDKN2B", "PTEN", "STK11"]

	level_12_var = [var for var in var_data if judge_var(var, [4,5], [4,5])]
	# 2025.02.21- 更新展示变异等级
	# 1. 有配对：体细胞I/II/肿瘤发生发展 + 胚系4/5 + 胚系3类但有用药
	# 2. 无配对：体细胞I/II/肿瘤发生发展 + 预测胚系4/5（归为I/II/肿瘤发生发展） + 预测胚系3类但有用药
	# ==> 即新增（确认/预测）胚系3类但有用药的变异即可
	germline_level3_regimen = [var for var in var_data if var["var_origin"] == "germline" and var["clinic_num_g"] == 3 and judgeRegimen(var)]
	if germline_level3_regimen:
		level_12_var.extend(germline_level3_regimen)
	# 2025.02.21-新增完成

	for var in level_12_var:
		# 仅展示扩增的基因
		if var["bio_category"] == "Cnv" and var["gene_symbol"] in cnv_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			io_result[var["gene_symbol"]].append(fdzs_mp_cnv_stran(var))
		# 仅展示融合的基因
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"] and set(re.split(",", var["gene_symbol"])) & set(["ALK"]):
			if "ALK" not in io_result.keys():
				io_result.setdefault("ALK", [])
			# 融合可能会出现exon相同断点不同的变异，报告中会重复，这边加个去重-2023.04.13
			if var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"融合" not in io_result["ALK"]:
				io_result["ALK"].append(var["five_prime_gene"]+":"+var["five_prime_cds"]+"-"+var["three_prime_gene"]+":"+var["three_prime_cds"]+"融合")
		# HD基因经确认同时展示snvindel和hd
		elif (var["bio_category"] == "Snvindel" or var["bio_category"] == "PHd") and var["gene_symbol"] in hd_gene_list:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["bio_category"] == "Snvindel":
				if var["hgvs_p"] != "p.?":
					io_result[var["gene_symbol"]].append(var["hgvs_p"])
				else:
					io_result[var["gene_symbol"]].append(var["hgvs_c"])
			elif var["bio_category"] == "PHd":
				if var["type"] == "HomoDel":
					if "纯合缺失" not in io_result[var["gene_symbol"]]:
						io_result[var["gene_symbol"]].append("纯合缺失")
				elif var["type"] == "HeteDel":
					if "杂合缺失" not in io_result[var["gene_symbol"]]:
						io_result[var["gene_symbol"]].append("杂合缺失")
				else:
					io_result[var["gene_symbol"]].append("未知变异类型！")
		# 其余基因展示Snvindel
		elif var["bio_category"] == "Snvindel" and var["gene_symbol"] in io_gene_P + io_gene_N:
			if var["gene_symbol"] not in io_result.keys():
				io_result.setdefault(var["gene_symbol"], [])
			if var["hgvs_p"] != "p.?":
				io_result[var["gene_symbol"]].append(var["hgvs_p"])
			else:
				io_result[var["gene_symbol"]].append(var["hgvs_c"])
	# summary展示
	both_cnv_list = ["CCND1","FGF3","FGF19"]
	io_p_list = ["{0} {1}".format(k, i) for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_P]
	io_n_list = ["{0} {1}".format(k, i) if not re.search("融合", i) else i for k,v in io_result.items() for i in v if k not in both_cnv_list and k in io_gene_N]

	# 处理CNV共突变
	if "CCND1" in io_result.keys() and "FGF3" in io_result.keys() and "FGF19" in io_result.keys():
		io_n_list.append("CCND1/FGF3/FGF19共扩增")

	return io_result, ", ".join(io_p_list), ", ".join(io_n_list)