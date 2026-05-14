#-*- coding:gbk -*-
import copy
from libs import listResultToDict

'''
Discription
	
	该脚本用来获取RNA表达结果，并根据模板需求进行分列处理。 

'''

def getRNA_exp(jsonDict):
	rna_exp_data = copy.deepcopy(jsonDict["rna_exp"])

	# 新增一个兼容-适配v4-2023.11.06
	# tpm为空时，v3返回“0.0”，v4返回变成了空值，反馈过说是要改报告脚本
	for i in rna_exp_data:
		if not i["tpm"]:
			i["tpm"] = "0.0"
	# 兼容完成-2023.11.06
	
	rna_exp = {}
	if rna_exp_data:
		# 临检通用版展示的是四列
		rna_exp["column_4"] = process_data(rna_exp_data, 4) 
		# 浙肿展示五列
		rna_exp["column_5"] = process_data(rna_exp_data, 5)
		# 浙肿需要展示GEP相关基因的表达
		rna_exp["gep"] = gep_tpm(rna_exp_data)
		# 厦门市一MP组织-综合，仅展示77个RNA EXP基因
		xmsy_rna_data = get_77gene(rna_exp_data)
		rna_exp["xmsy_77gene"] = process_data(xmsy_rna_data, 4)
		# 复旦中山展示七列-2024.05.09
		rna_exp["column_7"] = process_data(rna_exp_data, 7)

	return rna_exp

# num指的是模板中要按几列展示，基因+TPM为一列
def process_data(rna_exp_data, num):
	rna_exp = []
	for i in range(0, len(rna_exp_data)-num, num):
		tmp_dict = {}
		for j in range(1, num + 1):
			tmp_dict["gene"+str(j)] = rna_exp_data[i+j-1]["gene_symbol"]
			tmp_dict["tpm"+str(j)] = rna_exp_data[i+j-1]["tpm"]
		rna_exp.append(tmp_dict)
	
	rest_gene = len(rna_exp_data) % num
	rest_tmp_dict = {}
	for j in range(1, num+1):
		rest_tmp_dict["gene"+str(j)] = ""
		rest_tmp_dict["tpm"+str(j)] = ""

	rest_num = 1
	last_row_num = len(rna_exp_data)-rest_gene if rest_gene != 0 else len(rna_exp_data)-rest_gene-num
	for j in range(last_row_num, len(rna_exp_data)):
		rest_tmp_dict["gene"+str(rest_num)] = rna_exp_data[j]["gene_symbol"]
		rest_tmp_dict["tpm"+str(rest_num)] = rna_exp_data[j]["tpm"]
		rest_num += 1
	rna_exp.append(rest_tmp_dict)

	return rna_exp

# GEP 相关基因的表达量-适用浙肿Master
def gep_tpm(rna_exp_data):
	gene_list = ["CCL5","CD27","CD274","CD276","CD8A","CMKLR1","CXCL9","CXCR6","HLA-DQA1","HLA-DRB1","HLA-E","IDO1","LAG3","NKG7","PDCD1LG2","PSMB10","STAT1","TIGIT"]
	gep_exp = {}
	for i in rna_exp_data:
		if i["gene_symbol"] in gene_list:
			gep_exp.setdefault(i["gene_symbol"].replace("-",""), "")
			gep_exp[i["gene_symbol"].replace("-","")] = i["tpm"]

	return gep_exp

# 适用厦门市一MP组织-综合，仅展示77个RNA EXP基因-2024.04.18
def get_77gene(rna_exp_data):
	gene_list = ['ALK', 'APC', 'ARID1A', 'ATF1', 'ATM', 'ATRX', 'BAP1', 'BCOR', 'BRAF', 'BRCA1', \
				 'BRCA2', 'BRIP1', 'CDH1', 'CDK4', 'CDKN2A', 'CDKN2B', 'CHEK1', 'CHEK2', 'CIC', 'CTNNB1', \
				 'DAXX', 'EGFR', 'EPCAM', 'ERBB2', 'ESR1', 'EWSR1', 'FANCL', 'FGF19', 'FGFR1', 'FGFR2', \
				 'FGFR3', 'FH', 'FLCN', 'FOXO1', 'HRAS', 'IDH1', 'IDH2', 'KDR', 'KIT', 'KRAS', \
				 'MAML2', 'MDM2', 'MEN1', 'MET', 'MGMT', 'MLH1', 'MSH2', 'MSH6', 'MUTYH', 'MYB', \
				 'MYC', 'NCOA1', 'NF1', 'NF2', 'NRAS', 'NTRK1', 'NTRK2', 'PDGFRA', 'PIK3CA', 'PMS2', \
				 'POLD1', 'PTEN', 'RAD51C', 'RAD54L', 'RB1', 'RET', 'SDHA', 'SETD2', 'SMAD4', 'SMARCA4', \
				 'STK11', 'TERT', 'TP53', 'TSC1', 'TSC2', 'VEGFA', 'VHL']
	xmsy_rna_data = [i for i in rna_exp_data if i["gene_symbol"] in gene_list]
	return xmsy_rna_data