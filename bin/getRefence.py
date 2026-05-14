#-*- coding:gbk -*-
import xlrd
import os
import re
from libs.getInterRef import getRef_from_inter
'''
Discription
	
	该脚本用来获取固定参考文献和动态参考文献。 

'''
#BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#fixed_refer_path = os.path.join(BASE_DIR, "config/reference.xlsx")

# 格式化配置文件
def stran_xlrd(sheet_name, config):
	fixed_refer_path = os.path.join(config, "report_requirenment.xlsx")
	xls = xlrd.open_workbook(fixed_refer_path)
	refer_sheet = xls.sheet_by_name(sheet_name)
	key = refer_sheet.row_values(0)
	Data = []
	for num in range(1, refer_sheet.nrows):
		rows = refer_sheet.row_values(num)
		tmpdict = {}
		for i in range(len(key)):
			tmpdict[key[i]] = rows[i]
		Data.append(tmpdict)
	return Data	

# 固定文献格式处理
def stran_fixed_refer(i):
	# authors 取前三，用","连接，无空格，后跟et al.
	authors = i['authors'].split(',')
	if len(authors) > 3:
		i['authors'] = ','.join([authors[0], authors[1], authors[2], 'et al.'])
	# date 格式为(2013)
	if i['date']:
		i['date'] = '(' + re.search(r'\d{4}', str(i["date"]).replace(".0", "")).group(0) + ')'
	# title 最后加"."
	if i['title']:
		i['title'] = i['title'].strip('.') + '.'
	# pmid 加[]
	if i['PMID']:
		i['PMID'] = ''.join(['[', 'PMID:', str(int(i['PMID'])), ']'])
	result = " ".join([i["authors"], i["date"], i["title"], i["journal"], i["vol"], i["PMID"]])

	return result


# 固定参考文献
# 还需要胃癌分子分型、GEP、TME参数
def getfixed_refer(report_name, tumor_list, KNB, MSI, config):
	Data = stran_xlrd("reference-fixed", config)
	# 匹配
	fixed_refer = []
	refer_type_dict = {}
	for i in [refer for refer in Data if refer["docx_template"] == report_name]:
		refer_type_dict.setdefault(i["refer_type"], [])
		refer_type_dict[i["refer_type"]].append(i)
	# 1. fixed：固定参考文献直接展示
	if "fixed" in refer_type_dict.keys():
		fixed_refer += refer_type_dict["fixed"]
	# 2. tumor: 只匹配癌种
	if "tumor" in refer_type_dict.keys():
		for refer in refer_type_dict["tumor"]:
			if refer["tumor_or_type"] in tumor_list:
				fixed_refer.append(refer)
	# 3. tumor_and_result: 匹配癌种和结果，目前有4种
	if "tumor_and_result" in refer_type_dict.keys():
		for refer in refer_type_dict["tumor_and_result"]:
			# 3.1. 胃癌分子分型GA_result
			#    若为胃癌，且tumor_or_type为GA_type，根据分型选择对应参考文献（GA_EB、GA_MSI、GA_stable）
			# 固定文献，不用区分分型
#			if "胃癌" in tumor_list and refer["tumor_or_type"] == "GA_type":
#				pass

			# 3.2. KNB
			if refer["tumor_or_type"] == "KNB" and KNB:
				fixed_refer.append(refer)

			# 3.3. GEP
			if "肺癌" in tumor_list and refer["tumor_or_type"] == "GEP":
				fixed_refer.append(refer)

			# 3.4. TME
			if "肺癌" in tumor_list and refer["tumor_or_type"] == "TME":
				fixed_refer.append(refer)

	# 4. result ： 仅匹配结果，目前有1种
	if "result" in refer_type_dict.keys():
		for refer in refer_type_dict["result"]:
			# 4.1. MSI-H
			if refer["tumor_or_type"] == "MSI-H" and "var_id" in MSI.keys() and MSI["var_id"] and MSI["var_id"] == "MSI-H":
				fixed_refer.append(refer)

	fixed_refer_stran = [stran_fixed_refer(i) for i in fixed_refer]
	
	return fixed_refer_stran

# 动态参考文献
# 将动态文献进行拆分，报告按需使用
def getdynamic_refer(jsonDict, var_data, msi, hrd, var_brca):
	dynamic_refer = {}
	# 默认参考文献
	# 1. s_var12：获取基因介绍、变异解读、治疗、诊断、预后、遗传风险参考文献，包含I/II类
	dynamic_refer["s_var12"] = []
	for var in var_data["var_somatic"]["level_I"] + var_data["var_somatic"]["level_II"]:
		# 2026.03.02-融合若两个基因都在检测范围，需要都提取PMID
		if var["bio_category"] in ["Sv", "PSeqRnaSv"] and "," in var["gene_symbol"]:
			dynamic_refer["s_var12"].extend(getRef_from_inter(jsonDict, var["five_prime_gene_function"]))
			dynamic_refer["s_var12"].extend(getRef_from_inter(jsonDict, var["three_prime_gene_function"]))
		else:
			dynamic_refer["s_var12"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		# 2026.03.02-更新完成
		dynamic_refer["s_var12"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
		dynamic_refer["s_var12"].extend(var["evi_sum"]["refer_evi"])
	# 2. s_var_onco_nodrug 获取基因介绍、变异解读，包含肿瘤发生发展相关变异
	dynamic_refer["s_var_onco_nodrug"] = []
	for var in var_data["var_somatic"]["level_onco_nodrug"]:
		# 2026.03.02-融合若两个基因都在检测范围，需要都提取PMID
		if var["bio_category"] in ["Sv", "PSeqRnaSv"] and "," in var["gene_symbol"]:
			dynamic_refer["s_var_onco_nodrug"].extend(getRef_from_inter(jsonDict, var["five_prime_gene_function"]))
			dynamic_refer["s_var_onco_nodrug"].extend(getRef_from_inter(jsonDict, var["three_prime_gene_function"]))
		else:
			dynamic_refer["s_var_onco_nodrug"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		# 2026.03.02-更新完成
		dynamic_refer["s_var_onco_nodrug"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
	# 3. s_var3：获取基因介绍、变异解读，包含III类变异
	dynamic_refer["s_var3"] = []
	for var in var_data["var_somatic"]["level_III"]:
		# 2026.03.02-融合若两个基因都在检测范围，需要都提取PMID
		if var["bio_category"] in ["Sv", "PSeqRnaSv"] and "," in var["gene_symbol"]:
			dynamic_refer["s_var3"].extend(getRef_from_inter(jsonDict, var["five_prime_gene_function"]))
			dynamic_refer["s_var3"].extend(getRef_from_inter(jsonDict, var["three_prime_gene_function"]))
		else:
			dynamic_refer["s_var3"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		# 2026.03.02-更新完成
		dynamic_refer["s_var3"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
	# 4. g_var12: 获取基因介绍、变异解读、治疗、诊断、预后、遗传风险参考文献
	dynamic_refer["g_var45"] = []
	for var in var_data["var_germline"]["level_4"] + var_data["var_germline"]["level_5"]:
		dynamic_refer["g_var45"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		dynamic_refer["g_var45"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
		dynamic_refer["g_var45"].extend(var["evi_sum"]["refer_evi"])

	### MLPA 文献漏加了-del加到g_var45中，刘炜芬-2023.11.07
	for var in var_brca["mlpa"]["B1_Loss"] + var_brca["mlpa"]["B2_Loss"]:
		dynamic_refer["g_var45"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		dynamic_refer["g_var45"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
		if var["evi_sum"] and "refer_evi" in var["evi_sum"].keys() and var["evi_sum"]["refer_evi"]:
			dynamic_refer["g_var45"].extend(var["evi_sum"]["refer_evi"])
	### MLPA del添加完成-2023.11.07

	# 5. g_var3: 获取基因介绍、变异解读
	dynamic_refer["g_var3"] = []
	for var in var_data["var_germline"]["level_3"]:
		dynamic_refer["g_var3"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		dynamic_refer["g_var3"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))

	### MLPA 文献漏加了-dup加到g_var3中，刘炜芬-2023.11.07
	for var in var_brca["mlpa"]["B1_Gain"] + var_brca["mlpa"]["B2_Gain"]:
		dynamic_refer["g_var3"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		dynamic_refer["g_var3"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
	### MLPA del添加完成-2023.11.07

	# 5. KNB：获取KNB阴性时的治疗、诊断、预后、遗传风险参考文献
	dynamic_refer["knb"] = var_data["knb"]["evi_sum"]["refer_evi"] if var_data["knb"] else []
	# 6. EC：EC分型临床辅助治疗决策
	dynamic_refer["ec_type"] = var_data["ec_type"]["evi_sum"]["refer_evi"] if var_data["ec_type"] and "evi_sum" in var_data["ec_type"].keys() else []
	# 7. MSI：MSI-H时的治疗策略
	dynamic_refer["msi"] = msi["evi_sum"]["refer_evi"] if msi and msi["var_id"] == "MSI-H" else []
	# 8. HRD：治疗策略
	dynamic_refer["hrd"] = hrd["evi_sum"]["refer_evi"] if hrd and hrd["evi_sum"] else []
	# 2023.06.08-hrd添加在gss中的参考文献
	if "gss" in var_data.keys() and var_data["gss"] and "gss" in var_data["gss"].keys() and var_data["gss"]["gss"] and "evi_sum" in var_data["gss"]["gss"].keys()\
		and var_data["gss"]["gss"]["evi_sum"]:
		dynamic_refer["hrd"].extend(var_data["gss"]["gss"]["evi_sum"]["refer_evi"])

	# 9. 孙逸仙HRD C，仅展示I/II类变异解读
	dynamic_refer["s_var12_syx"] = []
	for var in var_data["var_somatic"]["level_I"] + var_data["var_somatic"]["level_II"]:
		dynamic_refer["s_var12_syx"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
	# 新增，孙逸仙III类变异解读-2023.05.16
	dynamic_refer["s_var3_syx"] = []
	for var in var_data["var_somatic"]["level_onco_nodrug"] + var_data["var_somatic"]["level_III"]:
		dynamic_refer["s_var3_syx"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
	# 新增，福建省立116 10 个基因III类基因介绍+变异解读-2023.05.24
	dynamic_refer["fjsl_116_lc10III"] = []
	for var in var_data["special"]["var_pan116_split"]["gene10_level_III"]:
		if var["bio_category"] == "Sv" and "," in var["gene_symbol"]:
			dynamic_refer["fjsl_116_lc10III"].extend(getRef_from_inter(jsonDict, var["five_prime_gene_function"]))
			dynamic_refer["fjsl_116_lc10III"].extend(getRef_from_inter(jsonDict, var["three_prime_gene_function"]))
		else:
			dynamic_refer["fjsl_116_lc10III"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		dynamic_refer["fjsl_116_lc10III"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))

	# 新增，厦门市一116，10基因和其他基因参考文献分开展示-2024.05.30
	dynamic_refer["xmsy_116_lc10"] = []
	for var in var_data["special"]["var_pan116_split"]["gene10_germline_45"] + var_data["special"]["var_pan116_split"]["gene10_level_I_II"]:
		if var["bio_category"] == "Sv" and "," in var["gene_symbol"]:
			dynamic_refer["xmsy_116_lc10"].extend(getRef_from_inter(jsonDict, var["five_prime_gene_function"]))
			dynamic_refer["xmsy_116_lc10"].extend(getRef_from_inter(jsonDict, var["three_prime_gene_function"]))
		else:
			dynamic_refer["xmsy_116_lc10"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		dynamic_refer["xmsy_116_lc10"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
		dynamic_refer["xmsy_116_lc10"].extend(var["evi_sum"]["refer_evi"])
		if var["var_origin"] == "germline":
			dynamic_refer["xmsy_116_lc10"].extend(var["evi_sum"]["refer_evi_risk"])

	dynamic_refer["xmsy_116_other"] = []
	for var in var_data["special"]["var_pan116_split"]["other_germline_45"] + var_data["special"]["var_pan116_split"]["other_level_I_II"]:
		if var["bio_category"] == "Sv" and "," in var["gene_symbol"]:
			dynamic_refer["xmsy_116_other"].extend(getRef_from_inter(jsonDict, var["five_prime_gene_function"]))
			dynamic_refer["xmsy_116_other"].extend(getRef_from_inter(jsonDict, var["three_prime_gene_function"]))
		else:
			dynamic_refer["xmsy_116_other"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		dynamic_refer["xmsy_116_other"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
		dynamic_refer["xmsy_116_other"].extend(var["evi_sum"]["refer_evi"])
		if var["var_origin"] == "germline":
			dynamic_refer["xmsy_116_other"].extend(var["evi_sum"]["refer_evi_risk"])
	# 2024.05.30-新增完成

	# 10. 分别获取胚系变异解读和基因简介中的参考文献， 4/5类变异还需提取遗传风险中的参考文献-2023.02.23
	dynamic_refer["g_var45_inter"] = []
	dynamic_refer["g_var45_gene"] = []
	dynamic_refer["g_var45_risk"] = []
	# 新增MLPA del-2023.11.07
	#for var in var_data["var_germline"]["level_5"] + var_data["var_germline"]["level_4"]:
	for var in var_data["var_germline"]["level_5"] + var_data["var_germline"]["level_4"] + var_brca["mlpa"]["B1_Loss"] + var_brca["mlpa"]["B2_Loss"]:
		dynamic_refer["g_var45_inter"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
		dynamic_refer["g_var45_gene"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		if var["evi_sum"] and "refer_evi_risk" in var["evi_sum"].keys() and var["evi_sum"]["refer_evi_risk"]:
			dynamic_refer["g_var45_risk"].extend(var["evi_sum"]["refer_evi_risk"])
	dynamic_refer["g_var3_inter"] = []
	dynamic_refer["g_var3_gene"] = []
	for var in var_data["var_germline"]["level_3"] + var_brca["mlpa"]["B1_Gain"] + var_brca["mlpa"]["B2_Gain"]:
		dynamic_refer["g_var3_inter"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
		dynamic_refer["g_var3_gene"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
	
	# cnv替换mlpa，新增一个参考文献-2024.09.11
	dynamic_refer["g_var45_inter_cnv"] = []
	dynamic_refer["g_var45_gene_cnv"] = []
	dynamic_refer["g_var45_risk_cnv"] = []
	for var in var_data["var_germline"]["level_5"] + var_data["var_germline"]["level_4"] + var_brca["gcnv"]["B1_Loss"] + var_brca["gcnv"]["B2_Loss"]:
		dynamic_refer["g_var45_inter_cnv"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
		dynamic_refer["g_var45_gene_cnv"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		dynamic_refer["g_var45_risk_cnv"].extend(var["evi_sum"]["refer_evi_risk"])
	dynamic_refer["g_var3_inter_cnv"] = []
	dynamic_refer["g_var3_gene_cnv"] = []
	for var in var_data["var_germline"]["level_3"] + var_brca["gcnv"]["B1_Gain"] + var_brca["gcnv"]["B2_Gain"]:
		dynamic_refer["g_var3_inter_cnv"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
		dynamic_refer["g_var3_gene_cnv"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
	# 新增完成-2024.09.11
	
	# 11. 获取I/II类变异的变异解读和基因介绍参考文献-2023.09.01
	dynamic_refer["s_var12_gene"] = []
	dynamic_refer["s_var12_inter"] = []
	for var in var_data["var_somatic"]["level_I"] + var_data["var_somatic"]["level_II"]:
		# 2026.03.02-融合若两个基因都在检测范围，需要都提取PMID
		if var["bio_category"] in ["Sv", "PSeqRnaSv"] and "," in var["gene_symbol"]:
			dynamic_refer["s_var12_gene"].extend(getRef_from_inter(jsonDict, var["five_prime_gene_function"]))
			dynamic_refer["s_var12_gene"].extend(getRef_from_inter(jsonDict, var["three_prime_gene_function"]))
		else:
			dynamic_refer["s_var12_gene"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		# 2026.03.02-更新完成
		dynamic_refer["s_var12_inter"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
	# 新增结束-2023.09.01
	# 将体细胞变异的基因介绍、变异解读和临床证据拆分开好了，模板中自由组合-2023.11.03
	dynamic_refer["s_var12_evi"] = []
	for var in var_data["var_somatic"]["level_I"] + var_data["var_somatic"]["level_II"]:
		dynamic_refer["s_var12_evi"].extend(var["evi_sum"]["refer_evi"])

	dynamic_refer["s_varonco_gene"] = []
	dynamic_refer["s_varonco_inter"] = []
	for var in var_data["var_somatic"]["level_onco_nodrug"]:
		# 2026.03.02-融合若两个基因都在检测范围，需要都提取PMID
		if var["bio_category"] in ["Sv", "PSeqRnaSv"] and "," in var["gene_symbol"]:
			dynamic_refer["s_varonco_gene"].extend(getRef_from_inter(jsonDict, var["five_prime_gene_function"]))
			dynamic_refer["s_varonco_gene"].extend(getRef_from_inter(jsonDict, var["three_prime_gene_function"]))
		else:
			dynamic_refer["s_varonco_gene"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		# 2026.03.02-更新完成
		dynamic_refer["s_varonco_inter"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
	
	dynamic_refer["s_var3_gene"] = []
	dynamic_refer["s_var3_inter"] = []
	for var in var_data["var_somatic"]["level_III"]:
		# 2026.03.02-融合若两个基因都在检测范围，需要都提取PMID
		if var["bio_category"] in ["Sv", "PSeqRnaSv"] and "," in var["gene_symbol"]:
			dynamic_refer["s_var3_gene"].extend(getRef_from_inter(jsonDict, var["five_prime_gene_function"]))
			dynamic_refer["s_var3_gene"].extend(getRef_from_inter(jsonDict, var["three_prime_gene_function"]))
		else:
			dynamic_refer["s_var3_gene"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		# 2026.03.02-更新完成
		dynamic_refer["s_var3_inter"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))

	# 胚系变异参考文献，改为三个内容，snvindel、mlpa和cnv，模板中按需组合-2024.12.05
	def germline_refer_class(var_list, class_type):
		refer_gene = []
		refer_inter = []
		refer_evi_sum = []
		refer_risk = []
		for var in var_list:
			refer_gene.extend(getRef_from_inter(jsonDict, var["gene_function"]))
			refer_inter.extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
			if var["evi_sum"] and "refer_evi" in var["evi_sum"].keys() and var["evi_sum"]["refer_evi"]:
				refer_evi_sum.extend(var["evi_sum"]["refer_evi"])
			if var["evi_sum"] and "refer_evi_risk" in var["evi_sum"].keys() and var["evi_sum"]["refer_evi_risk"]:
				refer_risk.extend(var["evi_sum"]["refer_evi_risk"])
		if class_type == "45":
			return refer_gene, refer_inter, refer_evi_sum, refer_risk
		else:
			return refer_gene, refer_inter
	# 1. snvindel-4/5类-分为基因介绍、变异解读、证据汇总（包含所有分类）、遗传风险四类
	g_snvindel_45 = var_data["var_germline"]["level_5"] + var_data["var_germline"]["level_4"]
	dynamic_refer["g_var45_gene_snvindel"], \
	dynamic_refer["g_var45_inter_snvindel"], \
	dynamic_refer["g_var45_evisum_snvindel"], \
	dynamic_refer["g_var45_risk_snvindel"] = germline_refer_class(g_snvindel_45, "45")
	# 2. snvindel-3类-分为基因介绍和变异解读两类
	g_snvindel_3 = var_data["var_germline"]["level_3"]
	dynamic_refer["g_var3_gene_snvindel"], \
	dynamic_refer["g_var3_inter_snvindel"] = germline_refer_class(g_snvindel_3, "3")
	# 3. MLPA-4/5类-分为基因介绍、变异解读、证据汇总（包含所有分类）、遗传风险四类
	mlpa_45 = var_brca["mlpa_v2"]["B1_mlpa_L5"] + var_brca["mlpa_v2"]["B2_mlpa_L5"] + var_brca["mlpa_v2"]["B1_mlpa_L4"] + var_brca["mlpa_v2"]["B2_mlpa_L4"]
	dynamic_refer["g_var45_gene_mlpa"], \
	dynamic_refer["g_var45_inter_mlpa"], \
	dynamic_refer["g_var45_evisum_mlpa"],\
	dynamic_refer["g_var45_risk_mlpa"] = germline_refer_class(mlpa_45, "45")
	# 4. mlpa-3类-分为基因介绍和变异解读两类
	mlpa_3 = var_brca["mlpa_v2"]["B1_mlpa_L3"] + var_brca["mlpa_v2"]["B2_mlpa_L3"]
	dynamic_refer["g_var3_gene_mlpa"], \
	dynamic_refer["g_var3_inter_mlpa"] = germline_refer_class(mlpa_3, "3")
	# 5. CNV-4/5类-分为基因介绍、变异解读、证据汇总（包含所有分类）、遗传风险四类
	g_cnv_45 = var_brca["gcnv_v2"]["B1_gcnv_L5"] + var_brca["gcnv_v2"]["B2_gcnv_L5"] + var_brca["gcnv_v2"]["B1_gcnv_L4"] + var_brca["gcnv_v2"]["B2_gcnv_L4"]
	dynamic_refer["g_var45_gene_cnv"], \
	dynamic_refer["g_var45_inter_cnv"], \
	dynamic_refer["g_var45_evisum_cnv"], \
	dynamic_refer["g_var45_risk_cnv"] = germline_refer_class(g_cnv_45, "45")
	# 6. CNV-3类-分为基因介绍和变异解读两类
	g_cnv_3 = var_brca["gcnv_v2"]["B1_gcnv_L3"] + var_brca["gcnv_v2"]["B2_gcnv_L3"]
	dynamic_refer["g_var3_gene_cnv"], \
	dynamic_refer["g_var3_inter_cnv"] = germline_refer_class(g_cnv_3, "3")
	# 2024.12.05-新增完成

	# 2025.03.12新增-gCNV参考文献
	dynamic_refer["germline_allcnv_var45_inter"] = []
	dynamic_refer["germline_allcnv_var45_gene"] = []
	dynamic_refer["germline_allcnv_var45_risk"] = []
	# 2026.01.28-新增evi_sum
	dynamic_refer["germline_allcnv_var45_evi"] = []
	# 2026.01.28-新增完成
	for var in var_data["germline_cnv"]["level_5"] + var_data["germline_cnv"]["level_4"]:
		dynamic_refer["germline_allcnv_var45_inter"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
		dynamic_refer["germline_allcnv_var45_gene"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		dynamic_refer["germline_allcnv_var45_risk"].extend(var["evi_sum"]["refer_evi_risk"])
		# 2026.01.28-新增evi_sum
		dynamic_refer["germline_allcnv_var45_evi"].extend(var["evi_sum"]["refer_evi"])
		# 2026.01.28-新增完成
	dynamic_refer["germline_allcnv_var3_inter"] = []
	dynamic_refer["germline_allcnv_var3_gene"] = []
	for var in var_data["germline_cnv"]["level_3"]:
		dynamic_refer["germline_allcnv_var3_inter"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
		dynamic_refer["germline_allcnv_var3_gene"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
	# 2025.03.12新增完成

	# 2025.12.24-同济MP非配对仅展示部分参考文献
	# 体系I/II（已有） + 胚系有用药位点
	dynamic_refer["tj_mp_nocontrol_germline_I_II_inter"] = []
	dynamic_refer["tj_mp_nocontrol_germline_I_II_gene"] = []
	dynamic_refer["tj_mp_nocontrol_germline_I_II_evi"] = []
	for var in var_data["var_germline"]["regimen_level_I"] + var_data["var_germline"]["regimen_level_II"]:
		dynamic_refer["tj_mp_nocontrol_germline_I_II_inter"].extend(getRef_from_inter(jsonDict, var["variant_interpret_cn"]))
		# 2026.03.02-融合若两个基因都在检测范围，需要都提取PMID
		if var["bio_category"] in ["Sv", "PSeqRnaSv"] and "," in var["gene_symbol"]:
			dynamic_refer["tj_mp_nocontrol_germline_I_II_gene"].extend(getRef_from_inter(jsonDict, var["five_prime_gene_function"]))
			dynamic_refer["tj_mp_nocontrol_germline_I_II_gene"].extend(getRef_from_inter(jsonDict, var["three_prime_gene_function"]))
		else:
			dynamic_refer["tj_mp_nocontrol_germline_I_II_gene"].extend(getRef_from_inter(jsonDict, var["gene_function"]))
		# 2026.03.02-更新完成
		dynamic_refer["tj_mp_nocontrol_germline_I_II_evi"].extend(var["evi_sum"]["refer_evi"])
	# 2025.12.24-新增完成
 
	return dynamic_refer