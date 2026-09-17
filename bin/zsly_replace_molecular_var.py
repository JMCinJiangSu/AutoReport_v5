#-*- coding:gbk -*-

'''
开发时间：2026.07.24
适用范围：中山大学附属第六医院 OncoPro（组织）

包含分子标志物：snvindel、cnv、sv、rna_sv、msi、hd、tmb、gss、mlpa、ec_type、knb。（其中cp200不包含rna_sv、tmb、gss、mlpa）
原始json文件->获取各个模块中变异var_code对应详细信息，字段参考之前发的示例->遍历各个模块中的molecular_var进行替换
使用修改后的json文件->正常格式转化->填充报告
'''


def var_code_stran(jsonDict):
    # oncopro用不上的先不加
    var_code_dict = {}
    # 1. snvindel
    for var in jsonDict["snvindel"]:
        var_code_dict[var["var_code"]]  = {
            "bio_category": "Snvindel",
            "gene_symbol": var["gene_symbol"],
            "hgvs_c": var["hgvs_c"],
            "hgvs_p": var["hgvs_p"],
            "gene_region": var["gene_region"],
            "freq": var["freq"]
        }
    # 2. cnv
    for var in jsonDict["cnv"]:
        var_code_dict[var["var_code"]]  = {
            "bio_category" : "Cnv",
            "gene_symbol" : var["gene_symbol"],
            "cn_mean" : var["cn_mean"]
        }
    # 3. sv
    for var in jsonDict["sv"]:
        var_code_dict[var["var_code"]] = {
            "bio_category" : "Sv",
            "gene_symbol" : var["gene_symbol"],
            "var_hgvs" : var["var_hgvs"],
            "copies" : var["copies"]
        }
    # 4. msi
    if jsonDict["msi"]:
        #print (jsonDict["msi"])
        var_code_dict[jsonDict["msi"][0]["var_code"]] = {
            "bio_category": "Special_markers",
            "biomarkers_name": jsonDict["msi"][0]["var_id"]
        }
    # 5. knb
    if jsonDict["knb"]:
        var_code_dict[jsonDict["knb"][0]["var_code"]] = {
            "bio_category": "Special_markers",
            "biomarkers_name": jsonDict["knb"][0]["var_id"]
        }
    # 6. hd
    if jsonDict["hd"]:
        for var in jsonDict["hd"]:
            var_code_dict[var["var_code"]] = {
                "bio_category": "PHd",
                "gene_symbol" : var["gene_symbol"],
                "region_exon" : var["region_exon"]
            }
    # 7. ec_type
    if jsonDict["ec_type"]:
        var_code_dict[jsonDict["ec_type"][0]["var_code"]] = {
            "bio_category": "Special_markers",
            "biomarkers_name": jsonDict["ec_type"][0]["var_id"]
        }
    return var_code_dict

def var_code_repalce(var, var_code_dict):
    #print (var)
    if var["evi_sum"]:
        for evi in var["evi_sum"]:
            molecular_var_list = []
            if evi["molecular_var"]:
                for var_code in evi["molecular_var"]:
                    molecular_var_list.append(var_code_dict.get(var_code, {"unknown" : var_code}))
            evi["molecular_var"] = molecular_var_list if molecular_var_list else []
    return var

def stran_json(jsonDict):
    var_code_list = var_code_stran(jsonDict)
    jsonDict["snvindel"] = [var_code_repalce(var, var_code_list) for var in jsonDict["snvindel"]]
    jsonDict["cnv"] = [var_code_repalce(var, var_code_list) for var in jsonDict["cnv"]]
    jsonDict["sv"] = [var_code_repalce(var, var_code_list) for var in jsonDict["sv"]]
    jsonDict["msi"] = [var_code_repalce(jsonDict["msi"][0], var_code_list)]
    jsonDict["knb"] = [var_code_repalce(jsonDict["knb"][0], var_code_list)] if jsonDict["knb"] else []
    # ec_type报告展示中没有归到I/II/III类变异中，这边就不用处理
    jsonDict["hd"] = [var_code_repalce(var, var_code_list) for var in jsonDict["hd"]]
    return jsonDict