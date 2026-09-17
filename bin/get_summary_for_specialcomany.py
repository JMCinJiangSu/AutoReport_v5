#-*- coding:gbk -*-
import re
from pypinyin import pinyin, Style
from itertools import chain

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

def ec_summary(ec_type):
	ec_dict = {
		"POLE-ultramutated type EC" : "POLE突变型（POLE mutation，POLE mut）",
		"MSI-H type EC" : "错配修复功能缺陷（Mismatch repair deficiency，MMRd）",
		"CNH type EC" : "TP53基因突变（p53 abnormality，p53 abn）",
		"CNL type EC" : "非特异性分子谱（Non-specific molecular profile，NSMP）"
	}
	return ec_dict.get(ec_type, ec_type)

def zjey_hrr_varsum(var_list):
	result = []
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			hgvs = var["hgvs_p"] if var["hgvs_p"] != "p.?" else var["hgvs_c"]
			freq_ = ""
			if var["var_origin"] == "germline":
				freq_ = "纯合" if float(var["freq"]) >= 0.85 else "杂合"
			else:
				freq_ = var["freq_str"]
			result.append("{0} {1}（{2}）".format(var["gene_symbol"], hgvs, freq_))
		elif var["type"] == "Loss":
			hete_or_homo = ""
			if "cnv_type" in var.keys() and var["cnv_type"]:
				hete_or_homo = "杂合" if var["cnv_type"] == "HeteDel" else "纯合" if var["cnv_type"] == "HomoDel" else ""
			if hete_or_homo:
				result.append("{0} {1} del（{2}）".format(var["gene_symbol"], var["value"], hete_or_homo))
			else:
				result.append("{0} {1} del".format(var["gene_symbol"], var["value"]))
		elif var["type"] == "Gain":
			result.append("{0} {1} dup".format(var["gene_symbol"], var["value"]))
	return "、".join(result)

def zjey_150_var_sum(var_list):
	result = []
	for var in var_list:
		if var["type"] == "Loss":
			result.append(var["gene_symbol"] + " " + var["value"] + " del")
		elif var["type"] == "Gain":
			result.append(var["gene_symbol"] + " " + var["value"] + " dup")
		else:
			if var["hgvs_p"] != "p.?":
				result.append(var["gene_symbol"]+" "+var["hgvs_c"]+" "+var["hgvs_p"])
			else:
				result.append(var["gene_symbol"]+" "+var["hgvs_c"])
	return ", ".join(result)

def zjfe_cp40_var_sum(var_list):
	result = []
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			if var["hgvs_p"] != "p.?":
				result.append("{0} {1}（{2}）".format(var["gene_symbol"], var["hgvs_p"], var["freq_str"]))
			else:
				result.append("{0} {1}（{2}）".format(var["gene_symbol"], var["hgvs_c"], var["freq_str"]))
		elif var["bio_category"] == "Cnv":
			result.append("{0} {1}（{2}）".format(var["gene_symbol"], "扩增", var["cn_mean"]))
		elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
			if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
				result.append("MET exon14 跳跃（{0} copies）".format(var["copies"]))
			else:
				result.append("{0}-{1} 融合（{2} copies）".format(var["five_prime_gene"], var["three_prime_gene"], var["copies"]))
	return ", ".join(result)

def var_s_summary(info):
	var_list = info[0]
	sample = info[1]
	def sum_var(vlist):
		v_result = []
		for var in vlist:
			if var["bio_category"] == "Snvindel":
				if var["hgvs_p"] != "p.?":
					v_result.append(var["gene_symbol"]+" "+var["hgvs_p"])
				else:
					v_result.append(var["gene_symbol"]+" "+var["hgvs_c"])
			elif var["bio_category"] == "Cnv":
				v_result.append(var["gene_symbol"]+" 扩增")
			elif var["bio_category"] in ["Sv", "PSeqRnaSv"]:
				if var["five_prime_gene"] == "MET" and var["three_prime_gene"] == "MET":
					v_result.append("MET exon14 跳跃")
				else:
					if var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合" not in v_result:
						v_result.append(var["five_prime_gene"]+"-"+var["three_prime_gene"]+" 融合")
			elif var["bio_category"] == "PHd":
				if var["type"] == "HomoDel":
					if var["gene_symbol"]+" 纯合缺失" not in v_result:
						v_result.append(var["gene_symbol"]+" 纯合缺失")
				elif var["type"] == "HeteDel":
					if var["gene_symbol"]+" 杂合缺失" not in v_result:
						v_result.append(var["gene_symbol"]+" 杂合缺失")
				else:
					v_result.append(var["gene_symbol"]+" 未知变异类型！")
		return ", ".join(v_result)
	c_var_all_num = len(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"] + \
					  var_list["var_somatic"]["level_onco_nodrug"] + var_list["var_somatic"]["level_III"])
	c_var_onco_drug_num = len(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"])
	c_var_onco_nodrug_num = len(var_list["var_somatic"]["level_onco_nodrug"])
	c_var_onco_drug_str = sum_var(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"])
	var_all_num = len(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"] + \
				 var_list["var_somatic"]["level_onco_nodrug"] + var_list["var_somatic"]["level_III"] +\
				 var_list["var_germline"]["level_5"] + var_list["var_germline"]["level_4"] + var_list["var_germline"]["regimen_noonco_var"])
	var_onco_drug_num = len(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"]+\
							var_list["var_germline"]["regimen_level_I"] + var_list["var_germline"]["regimen_level_II"])
	var_onco_nodrug_num = len(var_list["var_somatic"]["level_onco_nodrug"] + var_list["var_germline_nodrug"])
	var_onco_drug_str = sum_var(var_list["var_somatic"]["level_I"] + var_list["var_somatic"]["level_II"]+\
								  var_list["var_germline"]["regimen_level_I"] + var_list["var_germline"]["regimen_level_II"])
	if sample["control_sample_id"]:
		str_1 = "检出{0}个体细胞变异，其中具有临床意义的变异有{1}个，肿瘤发生发展相关变异有{2}个。".format(str(c_var_all_num), str(c_var_onco_drug_num), str(c_var_onco_nodrug_num))
		str_2 = "具有临床意义的变异有{0}。".format(c_var_onco_drug_str)
		if c_var_all_num != 0:
			if c_var_onco_drug_num != 0:
				return str_1 + str_2
			else:
				return str_1
		else:
			return "在检测范围内，未检出体细胞变异。"
	else:
		str_1 = "检出{0}个变异，其中具有临床意义的变异有{1}个，肿瘤发生发展相关变异有{2}个。".format(str(var_all_num), str(var_onco_drug_num), str(var_onco_nodrug_num))
		str_2 = "具有临床意义的变异有{0}。".format(var_onco_drug_str)
		if var_all_num != 0:
			if var_onco_drug_num != 0:
				return str_1 + str_2
			else:
				return str_1
		else:
			return "在检测范围内，未检出变异。"

def tme_summary(tme):
	tme_dict = {
	"IE/F" : "免疫富集/纤维化亚型(IE/F)",
	"IE" : "免疫富集/非纤维化亚型(IE)",
	"F" : "纤维化亚型(F)",
	"D" : "免疫荒漠型(D)"
	}
	return "TME分型为{0}".format(tme_dict.get(tme["tme_type"], tme["tme_type"]))

def io_summary_v4(info):
	io = info[0]
	sample = info[1]
	result = []
	if sample["prod_names"] in ["Master Panel（组织）", "Master Panel（血液）", "Master Panel V1（组织）"]:
		if io["hd_io_p_summary"]:
			result.append(io["hd_io_p_summary"]+"（疗效正相关）")
		if io["hd_io_n_summary"]:
			result.append(io["hd_io_n_summary"]+"（疗效负相关）")
	else:
		if io["io_p_summary_116"]:
			result.append(io["io_p_summary_116"]+"（疗效正相关）")
		if io["io_n_summary_116"]:
			result.append(io["io_n_summary_116"]+"（疗效负相关）")
	if not result:
		result = ["未检出相关变异"]
	return "；\n".join(result)+"。"

def hrd_summary_v2(info):
	gss = info[0]
	sample = info[1]
	result = []
	hrd_result = "HRD阳性" if gss["BRCA1"] or gss["BRCA2"] or float(gss["gss"]["gsscore"]) >= 45 else "HRD阴性"
	note = "（结果仅供参考）" if not gss["BRCA1"] and not gss["BRCA2"] and \
							(("var_auto_result" in gss["gss"].keys() and gss["gss"]["var_auto_result"] and gss["gss"]["var_auto_result"] == "F") or \
							 not sample["tumor_content"] or \
							 (sample["tumor_content"] and not is_number(sample["tumor_content_num"])) or \
							 (sample["tumor_content"] and is_number(sample["tumor_content_num"]) and float(sample["tumor_content_num"]) < 30)) \
							else ""
	result.append(hrd_result + note)
	if gss["summary"]:
		result.append("HRR通路相关基因突变：{0}".format(gss["summary"]))
	return "；\n".join(result)+"。"

def ga_summary(info):
	ga = info[0]
	msi = info[1]
	result = []
	if ga["ebv_type"]["ebv_type"] == "P" and msi["var_id"] == "MSI-H":
		if ga["ebv_sum"]:
			result.append("EB病毒感染型（EBV），{0}".format(ga["ebv_sum"]))
		else:
			result.append("EB病毒感染型（EBV）")
			result.append("微卫星不稳定型（MSI）")
	elif ga["ebv_type"]["ebv_type"] == "P":
		result.append("EB病毒感染型（EBV）")
	elif msi["var_id"] == "MSI-H":
		result.append("微卫星不稳定型（MSI）")
	else:
		if ga["gs_sum"]:
			result.append("基因组稳定型（GS）")
		if ga["cin_sum"]:
			result.append("染色体不稳定型（CIN）")
	if not result:
		result = ["-"]
	return "\n".join(result)

def var_g_summary(var_list):
	g_var_list = var_list["var_germline"]["level_5"] + var_list["var_germline"]["level_4"]
	result = []
	if g_var_list:
		for var in g_var_list:
			hgvs = var["hgvs_c"] if var["hgvs_p"] == "p.?" else var["hgvs_p"]
			clinic_g_cn = "致病性变异" if var["clinic_num_g"] == 5 else "疑似致病性变异"
			result.append("检出{0} {1}，为{2}".format(var["gene_symbol"], hgvs, clinic_g_cn))
	else:
		result = ["在检测范围内，未检出致病性/疑似致病性变异"]
	return "；".join(result)+"。"

def ZJFE_summary(data, judge_brca_cnv):
	germline_5 = data["var"]["var_germline"]["level_5"] + data["var_brca"]["mlpa_v2"]["B1_mlpa_L5"] + data["var_brca"]["mlpa_v2"]["B2_mlpa_L5"] if judge_brca_cnv == "mlpa" else \
		         data["var"]["var_germline"]["level_5"] + data["var_brca"]["gcnv_v2"]["B1_gcnv_L5"] + data["var_brca"]["gcnv_v2"]["B2_gcnv_L5"]
	germline_4 = data["var"]["var_germline"]["level_4"] + data["var_brca"]["mlpa_v2"]["B1_mlpa_L4"] + data["var_brca"]["mlpa_v2"]["B2_mlpa_L4"] if judge_brca_cnv == "mlpa" else \
		         data["var"]["var_germline"]["level_4"] + data["var_brca"]["gcnv_v2"]["B1_gcnv_L4"] + data["var_brca"]["gcnv_v2"]["B2_gcnv_L4"]
	germline_3 = data["var"]["var_germline"]["level_3"] + data["var_brca"]["mlpa_v2"]["B1_mlpa_L3"] + data["var_brca"]["mlpa_v2"]["B2_mlpa_L3"] if judge_brca_cnv == "mlpa" else \
		         data["var"]["var_germline"]["level_3"] + data["var_brca"]["gcnv_v2"]["B1_gcnv_L3"] + data["var_brca"]["gcnv_v2"]["B2_gcnv_L3"]
	somatic_I = data["var"]["var_somatic"]["level_I"]
	somatic_II = data["var"]["var_somatic"]["level_II"]
	somatic_onco_nodrug = data["var"]["var_somatic"]["level_onco_nodrug"]
	somatic_III = data["var"]["var_somatic"]["level_III"]
	gene150_germline_5 = data["var"]["var_germline"]["level_5"] + data["var"]["germline_cnv"]["level_5"]
	gene150_germline_4 = data["var"]["var_germline"]["level_4"] + data["var"]["germline_cnv"]["level_4"]
	gene150_germline_3 = data["var"]["var_germline"]["level_3"] + data["var"]["germline_cnv"]["level_3"]

	if data["sample"]["prod_names"] == "HRR（全血）":
		if germline_5 + germline_4 + germline_3:
			if germline_5:
				germline_var_sum = "检出{0}个变异，其中致病性变异有{1}个，疑似致病性变异有{2}个，意义不明确变异有{3}个。具有致病性变异有{4}。".format(
							len(germline_5 + germline_4 + germline_3), len(germline_5), len(germline_4), len(germline_3), zjey_hrr_varsum(germline_5))
			else:
				germline_var_sum = "检出{0}个变异，其中致病性变异有{1}个，疑似致病性变异有{2}个，意义不明确变异有{3}个。".format(
							len(germline_5 + germline_4 + germline_3), len(germline_5), len(germline_4), len(germline_3))
		else:
			germline_var_sum = "在检测范围内，未检测到基因变异。"
		result = {"germline_var_sum" : germline_var_sum}
			
	elif data["sample"]["prod_names"] == "HRR（组织）":
		if somatic_I + somatic_II + somatic_onco_nodrug + somatic_III:
			if somatic_I + somatic_II:
				unknow_origin_var_sum = "检出{0}个变异，其中具有临床意义的变异有{1}个，可能与肿瘤发生发展相关变异有{2}个。具有临床意义的变异有{3}。".format(
							len(somatic_I + somatic_II + somatic_onco_nodrug + somatic_III), len(somatic_I + somatic_II), len(somatic_onco_nodrug), zjey_hrr_varsum(somatic_I + somatic_II))
			else:
				unknow_origin_var_sum = "检出{0}个变异，其中具有临床意义的变异有{1}个，可能与肿瘤发生发展相关变异有{2}个。".format(
							len(somatic_I + somatic_II + somatic_onco_nodrug + somatic_III), len(somatic_I + somatic_II), len(somatic_onco_nodrug))
		else:
			unknow_origin_var_sum = "在检测范围内，未检测到基因变异。"
		result = {"unknow_origin_var_sum" : unknow_origin_var_sum}
			
	elif data["sample"]["prod_names"] == "HRR（组织 全血）":
		if somatic_I + somatic_II + somatic_onco_nodrug + somatic_III:
			if somatic_I + somatic_II:
				somatic_var_sum = "检出{0}个体细胞变异，其中具有临床意义的变异有{1}个，可能与肿瘤发生发展相关变异有{2}个。具有临床意义的变异有{3}。".format(
							len(somatic_I + somatic_II + somatic_onco_nodrug + somatic_III), len(somatic_I + somatic_II), len(somatic_onco_nodrug), zjey_hrr_varsum(somatic_I + somatic_II))
			else:
				somatic_var_sum = "检出{0}个体细胞变异，其中具有临床意义的变异有{1}个，可能与肿瘤发生发展相关变异有{2}个。".format(
							len(somatic_I + somatic_II + somatic_onco_nodrug + somatic_III), len(somatic_I + somatic_II), len(somatic_onco_nodrug))
		else:
			somatic_var_sum = "在检测范围内，未检测到体细胞变异。"
			
		if germline_5 + germline_4 + germline_3:
			if germline_5:
				germline_var_sum = "检出{0}个胚系变异，其中致病性变异有{1}个，疑似致病性变异有{2}个，意义不明确变异有{3}个。具有致病性变异有{4}。".format(
							len(germline_5 + germline_4 + germline_3), len(germline_5), len(germline_4), len(germline_3), zjey_hrr_varsum(germline_5))
			else:
				germline_var_sum = "检出{0}个胚系变异，其中致病性变异有{1}个，疑似致病性变异有{2}个，意义不明确变异有{3}个。".format(
							len(germline_5 + germline_4 + germline_3), len(germline_5), len(germline_4), len(germline_3))
		else:
			germline_var_sum = "在检测范围内，未检测到胚系变异。"
		result = {
			"somatic_var_sum" : somatic_var_sum,
			"germline_var_sum" : germline_var_sum
        }

	# 2026.08.28-150总结描述有调整		
	elif data["sample"]["prod_names"] == "遗传易感150基因":
		if gene150_germline_5 + gene150_germline_4 + gene150_germline_3:
			tmp_list = []
			if gene150_germline_5:
				#tmp_list.append("具有致病性的变异有{0}".format(zjey_150_var_sum(gene150_germline_5)))
				tmp_list.append("致病性变异为{0}".format(zjey_150_var_sum(gene150_germline_5)))
			if gene150_germline_4:
				#tmp_list.append("具有疑似致病性的变异有{0}".format(zjey_150_var_sum(gene150_germline_4)))
				tmp_list.append("疑似致病性变异为{0}".format(zjey_150_var_sum(gene150_germline_4)))
			if gene150_germline_5 + gene150_germline_4:
				#germline_var_sum = "检出{0}个变异 ，其中具有致病性的变异有{1}个，疑似致病性的变异有{2}个。{3}。".format(
				#			len(gene150_germline_5 + gene150_germline_4 + gene150_germline_3), len(gene150_germline_5), len(gene150_germline_4), "；".join(tmp_list))
				germline_var_sum = "检出{0}个变异 ，其中致病性变异有{1}个，疑似致病性变异有{2}个。{3}。".format(
							len(gene150_germline_5 + gene150_germline_4 + gene150_germline_3), len(gene150_germline_5), len(gene150_germline_4), "；".join(tmp_list))
			else:
				#germline_var_sum = "检出{0}个变异 ，其中具有致病性的变异有{1}个，疑似致病性的变异有{2}个。".format(
				#			len(gene150_germline_5 + gene150_germline_4 + gene150_germline_3), len(gene150_germline_5), len(gene150_germline_4))
				germline_var_sum = "检出{0}个变异 ，其中致病性变异有{1}个，疑似致病性变异有{2}个。".format(
							len(gene150_germline_5 + gene150_germline_4 + gene150_germline_3), len(gene150_germline_5), len(gene150_germline_4))
		else:
			germline_var_sum = "未检测到致病性、疑似致病性和意义不明确的胚系变异"
		result = {"germline_var_sum" : germline_var_sum}

	elif data["sample"]["prod_names"] == "Classic Panel":
		if somatic_I + somatic_II + somatic_onco_nodrug + somatic_III:
				if somatic_I + somatic_II:
					unknow_origin_var_sum = "检出{0}个变异，其中具有临床意义的变异有{1}个，肿瘤发生发展相关变异有{2}个。具有临床意义的变异有{3}。".format(
								len(somatic_I + somatic_II + somatic_onco_nodrug + somatic_III), len(somatic_I + somatic_II), len(somatic_onco_nodrug), zjfe_cp40_var_sum(somatic_I + somatic_II))
				else:
					unknow_origin_var_sum = "检出{0}个变异，其中具有临床意义的变异有{1}个，肿瘤发生发展相关变异有{2}个。".format(
								len(somatic_I + somatic_II + somatic_onco_nodrug + somatic_III), len(somatic_I + somatic_II), len(somatic_onco_nodrug))
		else:
			unknow_origin_var_sum = "在检测范围内，未检出基因变异。"
		result = {}
		result["unknow_origin_var_sum"] = unknow_origin_var_sum
		result["msi"] = "微卫星稳定（MSS）。" if data["msi"]["var_id"] == "MSS" else "微卫星不稳定（MSI-H）。"
		result["ec_type"] = ec_summary(data["var"]["ec_type"]["var_id"]) + "。" if "子宫内膜癌" in data["sample"]["tumor_list"] and "实体瘤" not in data["sample"]["tumor_names_cn"] else ""
		
	elif data["sample"]["prod_names"] == "Master Panel（组织）":
		result = {}
		result["unknow_origin_var_sum"] = var_s_summary([data["var"], data["sample"]]) if not data["sample"]["control_sample_id"] else ""
		result["somatic_var_sum"] = var_s_summary([data["var"], data["sample"]]) if data["sample"]["control_sample_id"] else ""
		result["germline_var_sum"] = var_g_summary(data["var"]) if data["sample"]["control_sample_id"] else ""
		result["msi"] = "微卫星稳定（MSS）。" if data["msi"]["var_id"] == "MSS" else "微卫星不稳定（MSI-H）。" if data["msi"]["var_id"] == "MSI-H" else "未返回MSI结果！"
		TMB_result = "低" if data["tmb"]["var_id"] == "TMB-L" else "高"
		result["tmb"] = "{0} Muts/Mb, 肿瘤突变负荷较{1}（{2}）。".format(data["tmb"]["TMB_value"], TMB_result, "TMB-H" if data["tmb"]["var_id"] == "TMB-H" else "TMB-L")
		result["gep_tme"] = "GEP分值为{0}分".format(data["gep"]["gep_score"]) + "；" + tme_summary(data["tme"]) + "。" if "肺癌" in data["sample"]["tumor_list"] and "实体瘤" not in data["sample"]["tumor_names_cn"] and data["sample"]["rna_sample_id"] else ""
		result["io"] = io_summary_v4([data["var"]["io"], data["sample"]])
		result["hrd"] = hrd_summary_v2([data["var"]["gss"], data["sample"]])
		result["ga_type"] = ga_summary([data["var"]["GA_type"], data["msi"]]) if "胃癌" in data["sample"]["tumor_list"] and "实体瘤" not in data["sample"]["tumor_names_cn"] and data["sample"]["rna_sample_id"] else ""
		result["ec_type"] = ec_summary(data["var"]["ec_type"]["var_id"]) if "子宫内膜癌" in data["sample"]["tumor_list"] and "实体瘤" not in data["sample"]["tumor_names_cn"] else ""

	elif data["sample"]["prod_names"] == "BRCA1/BRCA2（全血）":
		result = []
		if judge_brca_cnv == "mlpa":
			gbrca_var_list = data["var_brca"]["snv_s"]["B1_L5"] + data["var_brca"]["mlpa_v2"]["B1_mlpa_L5"] + \
				             data["var_brca"]["snv_s"]["B2_L5"] + data["var_brca"]["mlpa_v2"]["B2_mlpa_L5"] + \
					         data["var_brca"]["snv_s"]["B1_L4"] + data["var_brca"]["mlpa_v2"]["B1_mlpa_L4"] + \
							 data["var_brca"]["snv_s"]["B2_L4"] + data["var_brca"]["mlpa_v2"]["B2_mlpa_L4"] + \
							 data["var_brca"]["snv_s"]["B1_L3"] + data["var_brca"]["mlpa_v2"]["B1_mlpa_L3"] + \
						     data["var_brca"]["snv_s"]["B2_L3"] + data["var_brca"]["mlpa_v2"]["B2_mlpa_L3"]
		else:
			gbrca_var_list = data["var_brca"]["snv_s"]["B1_L5"] + data["var_brca"]["gcnv_v2"]["B1_gcnv_L5"] + \
		                     data["var_brca"]["snv_s"]["B2_L5"] + data["var_brca"]["gcnv_v2"]["B2_gcnv_L5"] + \
                             data["var_brca"]["snv_s"]["B1_L4"] + data["var_brca"]["gcnv_v2"]["B1_gcnv_L4"] + \
						     data["var_brca"]["snv_s"]["B2_L4"] + data["var_brca"]["gcnv_v2"]["B2_gcnv_L4"] + \
						     data["var_brca"]["snv_s"]["B1_L3"] + data["var_brca"]["gcnv_v2"]["B1_gcnv_L3"] + \
						     data["var_brca"]["snv_s"]["B2_L3"] + data["var_brca"]["gcnv_v2"]["B2_gcnv_L3"]
		for var in gbrca_var_list:
			if var["type"] == "Loss":
				var_info = var["value"] + " del"
				freq = "杂合" if var["cnv_type"] == "HeteDel" else "纯合" if var["cnv_type"] == "HomoDel" else "/"
			elif var["type"] == "Gain":
				var_info = var["value"] + " dup"
				freq = "/"
			else:
				if var["hgvs_p"] != "p.?":
					var_info = var["gene_region"] + " " + var["hgvs_c"] + " " + var["hgvs_p"]
				else:
					var_info = var["gene_region"] + " " + var["hgvs_c"]
				freq = "纯合" if float(var["freq"]) >= 0.85 else "杂合"
			result.append({
				"gene_symbol" : var["gene_symbol"],
				"var_info" : var_info,
				"var_origin" : "胚系",
				"freq" : freq,
				"var_level" : "致病性" if var["clinic_num_g"] == 5 else "疑似致病性" if var["clinic_num_g"] == 4 else "意义不明确"
            })
			
	elif data["sample"]["prod_names"] == "BRCA1/BRCA2（组织）":
		result = []
		tbrca_var_list = data["var_brca"]["snv_s"]["B1_L5"] + data["var_brca"]["snv_s"]["B2_L5"] + \
					     data["var_brca"]["snv_s"]["B1_L4"] + data["var_brca"]["snv_s"]["B2_L4"] + \
						 data["var_brca"]["snv_s"]["B1_L3"] + data["var_brca"]["snv_s"]["B2_L3"]
		for var in tbrca_var_list:
			regimen = []
			if var["hgvs_p"] != "p.?":
				var_info = var["gene_region"] + " " + var["hgvs_c"] + " " + var["hgvs_p"]
			else:
				var_info = var["gene_region"] + " " + var["hgvs_c"]
			if "Predictive" in var["evi_sum"]["evi_split"]:
				for i in var["evi_sum"]["evi_split"]["Predictive"]:
					regimen.append("{0}（{1}, {2}级）".format(i["regimen_name"], i["clinical_significance_cn"], i["evi_conclusion_simple"]))
			result.append({
				"gene_symbol" : var["gene_symbol"],
				"var_info" : var_info,
				"var_origin" : "待定",
				"freq" : var["freq_str"],
				"var_level" : "I类-强临床意义" if var["clinic_num_s"] == 5 else "II类-潜在临床意义" if var["clinic_num_s"] == 4 else "III类-临床意义不明",
				"regimen" : regimen
            })
	# 2026.02.25-新增BPTM Plus组织
	elif data["sample"]["prod_names"] == "BPTM Plus（组织）":
		if somatic_I + somatic_II + somatic_onco_nodrug + somatic_III:
				if somatic_I + somatic_II:
					unknow_origin_var_sum = "检出{0}个变异，其中具有临床意义的变异有{1}个，可能与肿瘤发生发展相关变异有{2}个。具有临床意义的变异有{3}。".format(
								len(somatic_I + somatic_II + somatic_onco_nodrug + somatic_III), len(somatic_I + somatic_II), len(somatic_onco_nodrug), zjfe_cp40_var_sum(somatic_I + somatic_II))
				else:
					unknow_origin_var_sum = "检出{0}个变异，其中具有临床意义的变异有{1}个，可能与肿瘤发生发展相关变异有{2}个。".format(
								len(somatic_I + somatic_II + somatic_onco_nodrug + somatic_III), len(somatic_I + somatic_II), len(somatic_onco_nodrug))
		else:
			unknow_origin_var_sum = "在检测范围内，未检出基因变异。"
		result = {}
		result["unknow_origin_var_sum"] = unknow_origin_var_sum
		result["msi"] = "微卫星稳定（MSS）" if data["msi"]["var_id"] == "MSS" else "微卫星不稳定（MSI-H）"
		result["ec_type"] = ec_summary(data["var"]["ec_type"]["var_id"]) if "子宫内膜癌" in data["sample"]["tumor_list"] else ""
	# 2026.02.25-新增完成
	# 2026.04.16-新增BPTM Plus全血
	elif data["sample"]["prod_names"] == "BPTM Plus（全血）":
		result = []
		# 2025.05.28-新增5个林奇基因CNV，取消MLPA
		var_list = data["var"]["var_germline"]["level_5"] + data["var"]["germline_cnv"]["level_5"] + \
				   data["var"]["var_germline"]["level_4"] + data["var"]["germline_cnv"]["level_4"] + \
				   data["var"]["var_germline"]["level_3"] + data["var"]["germline_cnv"]["level_3"]
		#if judge_brca_cnv == "mlpa":
		#	var_list = data["var"]["var_germline"]["level_5"] + data["var_brca"]["mlpa_v2"]["B1_mlpa_L5"] + data["var_brca"]["mlpa_v2"]["B2_mlpa_L5"] + \
		#			   data["var"]["var_germline"]["level_4"] + data["var_brca"]["mlpa_v2"]["B1_mlpa_L4"] + data["var_brca"]["mlpa_v2"]["B2_mlpa_L4"] + \
		#			   data["var"]["var_germline"]["level_3"] + data["var_brca"]["mlpa_v2"]["B1_mlpa_L3"] + data["var_brca"]["mlpa_v2"]["B2_mlpa_L3"]
		#else:
		#	var_list = data["var"]["var_germline"]["level_5"] + data["var_brca"]["gcnv_v2"]["B1_gcnv_L5"] + data["var_brca"]["gcnv_v2"]["B2_gcnv_L5"] + \
		#			   data["var"]["var_germline"]["level_4"] + data["var_brca"]["gcnv_v2"]["B1_gcnv_L4"] + data["var_brca"]["gcnv_v2"]["B2_gcnv_L4"] + \
		#			   data["var"]["var_germline"]["level_3"] + data["var_brca"]["gcnv_v2"]["B1_gcnv_L3"] + data["var_brca"]["gcnv_v2"]["B2_gcnv_L3"]
		# 2025.05.28-更新完成
		for var in var_list:
			if var["type"] == "Loss":
				var_info = var["value"] + " del"
				freq = "杂合" if "cnv_type" in var.keys() and var["cnv_type"] == "HeteDel" else "纯合" if "cnv_type" in var.keys() and var["cnv_type"] == "HomoDel" else "/"
			elif var["type"] == "Gain":
				var_info = var["value"] + " dup"
				freq = "/"
			else:
				if var["hgvs_p"] != "p.?":
					var_info = var["gene_region"] + " " + var["hgvs_c"] + " " + var["hgvs_p"] + "\n" + var["transcript_primary"]
				else:
					var_info = var["gene_region"] + " " + var["hgvs_c"] + "\n" + var["transcript_primary"]
				freq = "纯合" if float(var["freq"]) >= 0.85 else "杂合"
			result.append(
				{
					"gene_symbol" : var["gene_symbol"],
					"var_info" : var_info,
					"freq" : freq,
					"var_level" : "致病性" if var["clinic_num_g"] == 5 else "疑似致病性" if var["clinic_num_g"] == 4 else "意义不明确"
				}
			)
	# 2026.04.16-新增完成
	
	return result

# 2026.06.15：新增福建协和tLC10结构化数据
def FJXH_tLC10_summary(data):
	var_list = data["var"]["var_somatic"]["level_I"] + data["var"]["var_somatic"]["level_II"] + data["var"]["var_somatic"]["level_onco_nodrug"] + data["var"]["var_somatic"]["level_III"]
	result = []
	for var in var_list:
		if var["bio_category"] == "Snvindel":
			result.append({
				"var_name" : "{0}:{1}:{2}:({3})".format(var["transcript_primary"], \
														var["gene_region"].replace("exon", "Exon").replace("intron", "Intron"), \
														var["hgvs_c"], \
														var["hgvs_p_ZJZL"]),
				"var_category" : "Snvindel",
				"content" : var["variant_desc_cn"] + var["variant_interpret_cn"] + "该突变提示患者" + fjfy_regimen_sum(fjxh_evi_sort(var["evi_sum"]["evi_split"]["Predictive"])) +"。" if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"] else var["variant_desc_cn"] + var["variant_interpret_cn"]
			})
		elif var["bio_category"] == "Cnv":
			region_exon = var["region_exon"] if var["region_exon"] else ""
			result.append({
				"var_name" : var["gene_symbol"] + "_" + region_exon,
				"var_category" : "Cnv",
				"content" : var["variant_desc_cn"] + var["variant_interpret_cn"] + "该突变提示患者" + fjfy_regimen_sum(fjxh_evi_sort(var["evi_sum"]["evi_split"]["Predictive"])) +"。" if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"] else var["variant_desc_cn"] + var["variant_interpret_cn"]
			})
		elif var["bio_category"] == "Sv":
			result.append({
				"var_name" : var["var_name"] if "var_name" in var.keys() and var["var_name"] else "",
				"var_category" : "Sv",
				"content" : var["variant_desc_cn"] + var["variant_interpret_cn"] + "该突变提示患者" + fjfy_regimen_sum(fjxh_evi_sort(var["evi_sum"]["evi_split"]["Predictive"])) +"。" if "Predictive" in var["evi_sum"]["evi_split"].keys() and var["evi_sum"]["evi_split"]["Predictive"] else var["variant_desc_cn"] + var["variant_interpret_cn"]
			})

	return result

# 福建附一CP200、116对治疗方案进行总结-2025.11.28
# 该突变提示【可能对XXX、XXX（证据等级为A），XXX、XXX（证据等级为B）敏感；可能对XXX（证据等级为C）耐药。】
def fjfy_regimen_sum(regimen_list):
	regimen_dict = {}
	for regimen in regimen_list:
		if (regimen["evi_conclusion_simple"], regimen["clinical_significance_cn"]) not in regimen_dict.keys():
			regimen_dict.setdefault((regimen["evi_conclusion_simple"], regimen["clinical_significance_cn"]), [])
		regimen_dict[(regimen["evi_conclusion_simple"], regimen["clinical_significance_cn"])].append(regimen["regimen_name"])
	
	result_sense = []
	result_resis = []
	result = []
	for level in ["A", "B", "C", "D"]:
		if (level, "敏感") in regimen_dict.keys():
			result_sense.append("{0}（证据等级为{1}）".format("、".join(regimen_dict[(level, "敏感")]), level))
		if (level, "耐药") in regimen_dict.keys():
			result_resis.append("{0}（证据等级为{1}）".format("、".join(regimen_dict[(level, "耐药")]), level))
	if result_sense:
		result.append("可能对{0}敏感".format("，".join(result_sense)))
	if result_resis:
		result.append("可能对{0}耐药".format("，".join(result_resis)))

	return "；".join(result)

# 福建协和tBRCA、tLC10临床提示排序更新-2025.01.15
# 同等级NMPA获批证据优先展示，FDA其次，再展示其他（获批字段看regimen_refer_agency）
# 再按首字母排
def get_regimen_first_str(regimen_name):
    regimen_first_str = []
    for regimen in re.split("\+", regimen_name):
        if "\u4e00" <= regimen <= "\u9fff":
            _str = []
            for i in regimen:
                if "\u4e00" <= i <= "\u9fff":
                    _str.append("".join(chain.from_iterable(pinyin(i, Style.TONE3)))[0].upper())
                else:
                    _str.append(i)
            regimen_first_str.extend(_str)
        else:
            regimen_first_str.append(regimen[0].upper())
    return "".join(regimen_first_str)

def fjxh_evi_sort(evi_list):
	for evi in evi_list:
		regimen_refer_agency = evi["regimen_refer_agency"] if "regimen_refer_agency" in evi.keys() and evi["regimen_refer_agency"] else ""
		if set(re.split(",", regimen_refer_agency)) & set(["NMPA"]):
			evi["fjxh_agency"] = 0
		elif set(re.split(",", regimen_refer_agency)) & set(["FDA"]):
			evi["fjxh_agency"] = 1
		else:
			evi["fjxh_agency"] = 2
		evi["regimen_name_FJXH_str"] = get_regimen_first_str(evi["regimen_name"])
	evi_list = sorted(evi_list, key = lambda i:(i["evi_conclusion_simple"], \
											  	i["sense_rule"], \
												i["fjxh_agency"], \
												i["regimen_name_FJXH_str"]))
	return evi_list

# 2026.08.07-福建肿瘤gBRCA结构化数据
def fjzl_gbrca_var_sum(data):
	var_list = data["var_brca"]["snv_s"]["B1_L5"] + data["var_brca"]["snv_s"]["B2_L5"] + \
			   data["var_brca"]["gcnv_v2"]["B1_gcnv_L5"] + data["var_brca"]["gcnv_v2"]["B2_gcnv_L5"] + \
			   data["var_brca"]["snv_s"]["B1_L4"] + data["var_brca"]["snv_s"]["B2_L4"] + \
			   data["var_brca"]["gcnv_v2"]["B1_gcnv_L4"] + data["var_brca"]["gcnv_v2"]["B2_gcnv_L4"]
	result = {
		"var" : [],
		"inter" : ""
	}
	for var in var_list:
		var_info = "{0}基因{1}号外显子发生大片段缺失变异".format(var["gene_symbol"], var["value"].replace("exon", "")) if var["type"] == "Loss" else \
			       "{0}基因{1}号外显子发生大片段重复变异".format(var["gene_symbol"], var["value"].replace("exon", "")) if var["type"] == "Gain" else \
				   "{0}基因{1}".format(var["gene_symbol"], var["varInter_FJZL"])
		var_inter_sum = "该患者外周血基因组DNA样品检测出{0}，{1}该位点证据等级为({2})，判读为{3}，预测其对铂类和PARP抑制剂的反应性可能增高。".format(
						var_info,
						var["variant_interpret_cn"],
						var["evidence_categorys"].replace(";", "+") if "evidence_categorys" in var.keys() and var["evidence_categorys"] else "NF",
						"致病性变异" if var["clinic_num_g"] == 5 else "疑似致病性变异"
					)
		result["var"].append({
			"var_name" : var["var_name"] if "var_name" in var.keys() and var["var_name"] else "",
			"interp" : var_inter_sum
		})
	var_sum = []
	for i in ["1", "2"]:
		L5_var = data["var_brca"]["snv_s"]["B" + i + "_L5"] + data["var_brca"]["gcnv_v2"]["B" + i + "_gcnv_L5"]
		L4_var = data["var_brca"]["snv_s"]["B" + i + "_L4"] +data["var_brca"]["gcnv_v2"]["B" + i + "_gcnv_L4"]
		if L5_var and L4_var:
			var_sum.append("BRCA" + i + "基因致病性和疑似致病性变异")
		elif L5_var:
			var_sum.append("BRCA" + i + "基因致病性变异")
		elif L4_var:
			var_sum.append("BRCA" + i + "基因疑似致病性变异")

	result["inter"] = "该患者外周血基因组DNA样品检出{0}，建议对患者进行遗传性乳腺癌/卵巢癌相关遗传咨询，并在适当情况下对其可能存在风险的直系亲属行该位点的胚系突变检测，至我院体检中心或专科医生处进行遗传咨询。".format("、".join(var_sum)) if var_sum else ""
	return result