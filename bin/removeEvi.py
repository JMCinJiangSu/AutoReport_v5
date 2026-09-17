#-*- coding:gbk -*-
import re

# 2025.08.19
# 后续有需要删除证据的特殊需求写在这里

def cqxn_remove_PIK3CA_evi(jsonDict):
    # 条件如下 #
    cqxn_prod_list = ["Pan116（组织）"]
    report_module_type = "hospital"
    company = "重庆西南医院"
    # 判断PIK3CA变异是否获批
    # 仅snvindel 11个位点获批，其他snvindel或cnv、sv都删除证据
    def judge_pik3ce_appr(var):
        result = False
        appr_list = ["C420R", "E542K", "E545A", "c.1635G>T", "E545G", "E545K", \
                     "Q546E", "Q546R", "H1047L", "H1047R", "H1047Y"]
        if var["bio_category"] == "Snvindel":
            if var["hgvs_c"] in appr_list:
                result = True
            elif var["hgvs_p"].replace("p.", "").replace("(", "").replace(")", "") in appr_list:
                result = True
        return result
    
    if jsonDict["sample_info"]["report_module_type"] == report_module_type and \
       jsonDict["sample_info"]["prod_names"] in cqxn_prod_list and \
       jsonDict["sample_info"]["company"] == company:
        for var in jsonDict["snvindel"]:
            if var["gene_symbol"] == "PIK3CA" and not judge_pik3ce_appr(var):
                var["evi_sum"] = []
        for var in jsonDict["cnv"]:
            if var["gene_symbol"] == "PIK3CA":
                var["evi_sum"] = []
        for var in jsonDict["sv"]:
            if "PIK3CA" in re.split(",", var["gene_symbol"]):
                var["evi_sum"] = []
    return jsonDict

# 西安交大一
# 10、116/18/21、mp、cp删除KRAS、NRAS中的呋喹替尼、瑞戈替尼、贝伐珠单抗（引用机构需要时CSCO）
def xajdy_remove_ras_evi(jsonDict):
    # 条件如下 #
    xajdy_prod_list = ["10基因（组织）", "10基因（血液）", "Pan116（组织）", "Pan116（血液）", "TC21（组织）", \
                      "TC21（血液）", "GA18（组织）", "GA18（血液）", "Master Panel（组织）", "Classic Panel"]
    report_module_type = "hospital"
    company = ["西安交通大学第一附属医院", "西安交通大学医学院第一附属医院"] 
    if jsonDict["sample_info"]["report_module_type"] == report_module_type and \
       jsonDict["sample_info"]["prod_names"] in xajdy_prod_list and \
       jsonDict["sample_info"]["company"] in company:
        for var in jsonDict["snvindel"]:
            # 2026.02.10-BRAF V600E也需要删除证据
            if var["gene_symbol"] in ["KRAS", "NRAS"] or var["gene_symbol"] == "BRAF" and var["hgvs_p"] in ["p.(V600E)", "p.V600E"]:
                var["evi_sum"] = [evi for evi in var["evi_sum"] if not ("regimen_name" in evi.keys() and evi["regimen_name"] and evi["regimen_name"] in ["呋喹替尼", "瑞戈非尼", "贝伐珠单抗"] and evi["refer_agency"] == "CSCO")]
    return jsonDict

# 2026.04.01-重庆西南删除CDKN2A Inactivating Mutation中的奥希替尼证据
# 只检测snvindel
def cqxn_remove_CDKN2A_evi(jsonDict):
    # 条件如下 #
    cqxn_prod_list = ["Pan116（组织）"]
    report_module_type = "hospital"
    company = "重庆西南医院"
    if jsonDict["sample_info"]["report_module_type"] == report_module_type and \
       jsonDict["sample_info"]["prod_names"] in cqxn_prod_list and \
       jsonDict["sample_info"]["company"] == company:
        for var in jsonDict["snvindel"]:
            if var["gene_symbol"] == "CDKN2A" and "var_category_names" in var.keys() and var["var_category_names"] and "CDKN2A Inactivating Mutation" in var["var_category_names"]:
                var["evi_sum"] = [evi for evi in var["evi_sum"] if not ("regimen_name" in evi.keys() and evi["regimen_name"] and evi["regimen_name"] == "奥希替尼")]
    return jsonDict

def cqxn_remove_KRAS_prognostic_evi(jsonDict):
    # 条件如下 #
    cqxn_prod_list = ["Pan116（组织）"]
    report_module_type = "hospital"
    company = "重庆西南医院"
    # 肺癌KRAS snvindel KRAS Activating Mutation 删除预后证据    
    if jsonDict["sample_info"]["report_module_type"] == report_module_type and \
       jsonDict["sample_info"]["prod_names"] in cqxn_prod_list and \
       jsonDict["sample_info"]["company"] == company:
        for var in jsonDict["snvindel"]:
            if var["gene_symbol"] == "KRAS" and "var_category_names" in var.keys() and var["var_category_names"] and "KRAS Activating Mutation" in var["var_category_names"] and \
                "肺癌" in jsonDict["sample_info"]["tumor_list"]:
                var["evi_sum"] = [evi for evi in var["evi_sum"] if evi["evidence_type"] != "Prognostic"]
    return jsonDict

# 中山市人民CP200删除结直肠癌PIK3CA PIK3R1 PTEN中的阿司匹林
def zsrm_cp200_remove_aspirin(jsonDict):
    # 条件如下 #
    zsrm_prod_list = ["OncoPro（组织）", "Classic Panel 200（组织）"]
    report_module_type = "hospital"
    company = "中山市人民医院"
    # 肠癌 PIK3CA Activating Mutation、PIK3R1 Inactivating Mutation、PTEN Inactivating Mutation删除阿司匹林(PTEN 0.1.4流程有HD)
    if jsonDict["sample_info"]["report_module_type"] == report_module_type and \
       jsonDict["sample_info"]["prod_names"] in zsrm_prod_list and \
       jsonDict["sample_info"]["company"] == company:
        for var in jsonDict["snvindel"] + jsonDict["hd"]:
            if "结直肠癌" in jsonDict["sample_info"]["tumor_list"] and \
                "var_category_names" in var.keys() and var["var_category_names"] and (\
                (var["gene_symbol"] == "PIK3CA" and "PIK3CA Activating Mutation" in var["var_category_names"]) or \
                (var["gene_symbol"] == "PIK3R1" and "PIK3R1 Inactivating Mutation" in var["var_category_names"]) or \
                (var["gene_symbol"] == "PTEN" and "PTEN Inactivating Mutation" in var["var_category_names"])
                ):
                var["evi_sum"] = [evi for evi in var["evi_sum"] if not ("regimen_name" in evi.keys() and evi["regimen_name"] and evi["regimen_name"] == "阿司匹林")]
    return jsonDict

# 广东医科附属CP200/116/HRR删除PIK3CA PIK3R1 PTEN中的阿司匹林（不限癌种）
def gdykfs_remove_aspirin(jsonDict):
    # 条件如下 #
    gdykfs_prod_list = ["OncoPro（组织）", "Classic Panel 200（组织）", "Pan116（血液）", "Pan116（组织）", "HRR（全血）", "HRR（组织）", "HRR（组织 全血）"]
    report_module_type = "hospital"
    company = "广东医科大学附属医院"
    # PIK3CA Activating Mutation、PIK3R1 Inactivating Mutation、PTEN Inactivating Mutation删除阿司匹林(PTEN 0.1.4流程有HD)
    if jsonDict["sample_info"]["report_module_type"] == report_module_type and \
       jsonDict["sample_info"]["prod_names"] in gdykfs_prod_list and \
       jsonDict["sample_info"]["company"] == company:
        for var in jsonDict["snvindel"] + jsonDict["hd"]:
            if "var_category_names" in var.keys() and var["var_category_names"] and (\
                (var["gene_symbol"] == "PIK3CA" and "PIK3CA Activating Mutation" in var["var_category_names"]) or \
                (var["gene_symbol"] == "PIK3R1" and "PIK3R1 Inactivating Mutation" in var["var_category_names"]) or \
                (var["gene_symbol"] == "PTEN" and "PTEN Inactivating Mutation" in var["var_category_names"])
                ):
                var["evi_sum"] = [evi for evi in var["evi_sum"] if not ("regimen_name" in evi.keys() and evi["regimen_name"] and evi["regimen_name"] == "阿司匹林")]
    return jsonDict

# 2026.06.30-武汉同济MP，不存在EGFR敏感突变时，KRAS不展示EGFR-TKIs耐药证据
def tj_mp_remove_KRAS_egfrtiks(jsonDict):
    # 条件如下 #
    tj_mp_prod_list = ["Master Panel（组织）"]
    report_module_type = "hospital"
    company = "华中科技大学同济医学院附属同济医院"
    # 判定是否存在EGFR敏感突变
    judge_egfr_var = False
    for var in jsonDict["snvindel"]:
        if ("var_category_names" in var.keys() and var["var_category_names"] and "EGFR Exon19 del" in var["var_category_names"]) or \
            var["hgvs_p"] in ["p.L858R", "p.L861Q", "p.S768I"] or \
            "G719" in var["hgvs_p"]:
            judge_egfr_var = True
            break

    if jsonDict["sample_info"]["report_module_type"] == report_module_type and \
       jsonDict["sample_info"]["prod_names"] in tj_mp_prod_list and \
       jsonDict["sample_info"]["company"] == company and not judge_egfr_var:
        for var in jsonDict["snvindel"]:
            if var["gene_symbol"] == "KRAS":
                var["evi_sum"] = [evi for evi in var["evi_sum"] if not ("regimen_name" in evi.keys() and evi["regimen_name"] and evi["regimen_name"] == "EGFR-TKIs")]
    return jsonDict

# 2026.07.03-西安交大一MP不展示MTAP HD变异
def xajdy_mp_remove_MTAP_HD(jsonDict):
    # 条件如下 #
    xajdy_mp_prod_list = ["Master Panel（组织）"]
    report_module_type = "hospital"
    company = "西安交通大学第一附属医院"
    # MTAP HD不展示
    if jsonDict["sample_info"]["report_module_type"] == report_module_type and \
       jsonDict["sample_info"]["prod_names"] in xajdy_mp_prod_list and \
       jsonDict["sample_info"]["company"] == company:
        jsonDict["hd"] = [var for var in jsonDict["hd"] if var["gene_symbol"] != "MTAP"]
    return jsonDict

# 2026.07.03-武汉同济MP
# 1. 用药、诊断、预后最高等级A/B时，删除预后、诊断和用药D证据
# 2. 用药、诊断、预后最高等级C/D时，删除预后、诊断证据
# 以上包含snvindel、cnv、sv、rna_sv、hd，其他分子标志物先不考虑
def whtj_mp_remove_evi(jsonDict):
     # 条件如下 #
    tj_mp_prod_list = ["Master Panel（组织）"]
    report_module_type = "hospital"
    company = "华中科技大学同济医学院附属同济医院"
    # I类过滤掉D级用药（辅助诊断和预后后面也是要过滤掉的，所以过滤的时候D级诊断/预后也过滤掉没关系）
    def filter_d(raw_evi_sum):
        level_list = [evi["evi_conclusion"][0] for evi in raw_evi_sum if evi["evidence_type"] in ["Predictive", "Prognostic", "Diagnostic"]]
        if set(["A", "B"]) & set(level_list):
            return [i for i in raw_evi_sum if i["evi_conclusion"][0] in ["A", "B", "C"]]
        else:
            return raw_evi_sum
    # 过滤掉预后和诊断
    def filter_dia_pro(raw_evi_sum):
        return [i for i in raw_evi_sum if i["evidence_type"] == "Predictive"]
    if jsonDict["sample_info"]["report_module_type"] == report_module_type and \
       jsonDict["sample_info"]["prod_names"] in tj_mp_prod_list and \
       jsonDict["sample_info"]["company"] == company:
        for var in jsonDict["snvindel"]:
            var["evi_sum"] = filter_dia_pro(filter_d(var["evi_sum"]))
        for var in jsonDict["cnv"]:
            var["evi_sum"] = filter_dia_pro(filter_d(var["evi_sum"]))
        for var in jsonDict["sv"]:
            var["evi_sum"] = filter_dia_pro(filter_d(var["evi_sum"]))
        for var in jsonDict["rna_sv"]:
            var["evi_sum"] = filter_dia_pro(filter_d(var["evi_sum"]))
        for var in jsonDict["hd"]:
            var["evi_sum"] = filter_dia_pro(filter_d(var["evi_sum"]))
    return jsonDict

# 2026.09.03-中山人民CP200、tHRR肺癌时删除KRAS预后证据
def zsrm_remove_KRAS_prognostic_evi(jsonDict):
    # 条件如下 #
    zsrm_prod_list = ["OncoPro（组织）", "Classic Panel 200（组织）", "HRR（组织）"]
    report_module_type = "hospital"
    company = "中山市人民医院"
    # 肺癌KRAS snvindel KRAS 删除预后证据    
    if jsonDict["sample_info"]["report_module_type"] == report_module_type and \
       jsonDict["sample_info"]["prod_names"] in zsrm_prod_list and \
       jsonDict["sample_info"]["company"] == company:
        for var in jsonDict["snvindel"]:
            if var["gene_symbol"] == "KRAS" and "肺癌" in jsonDict["sample_info"]["tumor_list"]:
                var["evi_sum"] = [evi for evi in var["evi_sum"] if evi["evidence_type"] != "Prognostic"]
    return jsonDict