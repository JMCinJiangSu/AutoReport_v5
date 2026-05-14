#-*- coding:gbk -*-

def remove_KRAS_evi(jsonDict):
    # 条件如下 #
    xajdy_prod_list = ["Pan116（组织）", "Pan116（血液）", "Master Panel（组织）", "Classic Panel", "10基因（组织）", "10基因（血液）","OncoPro（组织）", "Classic Panel 200（组织）"]
    report_module_type = "hospital"
    company = "西安交通大学第一附属医院"
    tumor = "肺癌"
    if jsonDict["sample_info"]["report_module_type"] == report_module_type and \
       jsonDict["sample_info"]["prod_names"] in xajdy_prod_list and \
       jsonDict["sample_info"]["company"] == company and \
       tumor in jsonDict["sample_info"]["tumor_list"]:
        for var in jsonDict["snvindel"]:       
            if var["gene_symbol"] == "KRAS" and var["evi_sum"]:
                A_prognostic_evi_sum = []
                evi_sum = []
                for evi in var["evi_sum"]:                
                    # KRAS EGFR-TKIs证据改为C级
                    if evi["regimen_name"] == "EGFR-TKIs" and "Resistant" in evi["clinical_significance"]:
                        print ("符合西安交大一相关条件，KRAS EGFR-TKIs证据改为C级，暂定C3！")
                        evi["evi_conclusion"] = "C3"
                        evi_sum.append(evi)
                    # 删除KRAS A级预后证据
                    elif evi["evidence_type"] == "Prognostic" and "A" in evi["evi_conclusion"]:
                        print ("符合西安交大一相关条件，KRAS A级预后证据不展示，不参与评级，仅在变异注释后面单独展示！")
                        if evi["evi_interpretation"] not in A_prognostic_evi_sum:
                            A_prognostic_evi_sum.append(evi["evi_interpretation"])
                    else:
                        evi_sum.append(evi)
                var["evi_sum"] = evi_sum
                # KRAS A级预后证据描述需要在变异注释后面单独展示，不参与证据评级
                var["xajdy_kras_a_prognostic"] = "".join(A_prognostic_evi_sum)
    return jsonDict