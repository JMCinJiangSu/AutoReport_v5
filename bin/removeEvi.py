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