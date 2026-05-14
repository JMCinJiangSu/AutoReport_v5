#-*- coding:gbk -*-

def remove_TP53_evi(jsonDict):
    # 条件如下 #
    xajdy_prod_list = ["Pan116（组织）", "Pan116（血液）", "Master Panel（组织）", "Classic Panel"]
    report_module_type = "hospital"
    company = "西安交通大学第一附属医院"
    tumor = "肺癌"
    if jsonDict["sample_info"]["report_module_type"] == report_module_type and \
       jsonDict["sample_info"]["prod_names"] in xajdy_prod_list and \
       jsonDict["sample_info"]["company"] == company and \
       tumor in jsonDict["sample_info"]["tumor_list"]:
        # TP53并且存在evi_sum，删除EGFR-TKIs证据
        for var in jsonDict["snvindel"]:
            if var["gene_symbol"] == "TP53" and var["evi_sum"]:
                for evi in var["evi_sum"]:
                    if evi["regimen_name"] == "EGFR-TKIs":
                        print ("符合西安交大一相关条件，删除TP53 EGFR-TKIs证据！")
                        var["evi_sum"].remove(evi) 
    return jsonDict

# 2025.06.12-重庆西南116 TP53不展示 迈华替尼
def cqzn_remove_TP53_evi(jsonDict):
    # 条件如下 #
    prod_list = ["Pan116（组织）"]
    report_module_type = "hospital"
    company = "重庆西南医院"
    if jsonDict["sample_info"]["report_module_type"] == report_module_type and \
       jsonDict["sample_info"]["prod_names"] in prod_list and \
       jsonDict["sample_info"]["company"] == company:
        # TP53存在迈华替尼 删除
        for var in jsonDict["snvindel"]:
            if var["gene_symbol"] == "TP53" and var["evi_sum"]:
                for evi in var["evi_sum"]:
                    # 2025.12.03-迈华替尼更新为美凡厄替尼
                    #if evi["regimen_name"] == "迈华替尼":
                    #    print ("符合重庆西南相关条件，删除TP53 迈华替尼证据！")
                    #    var["evi_sum"].remove(evi)
                    if evi["regimen_name"] == "美凡厄替尼":
                        print ("符合重庆西南相关条件，删除TP53 美凡厄替尼证据！")
                        var["evi_sum"].remove(evi) 
    return jsonDict