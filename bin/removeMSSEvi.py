#-*- coding:gbk -*-

def remove_MSS_evi(jsonDict):
    # 删除msi中MSS的证据
    if jsonDict["msi"]:
        if type(jsonDict["msi"]).__name__ == "dict":
            if "var_id" in jsonDict["msi"].keys() and jsonDict["msi"]["var_id"] and jsonDict["msi"]["var_id"] == "MSS":
                if "evi_sum" in jsonDict["msi"].keys():
                    jsonDict["msi"]["evi_sum"] = []
        else:
            if "var_id" in jsonDict["msi"][0].keys() and jsonDict["msi"][0]["var_id"] and jsonDict["msi"][0]["var_id"] == "MSS":
                if "evi_sum" in jsonDict["msi"][0].keys():
                    jsonDict["msi"][0]["evi_sum"] = []

    # 删除治疗方案介绍中的MSS therapeutic_regimen
    mss_biomarker = {"biomarker_type" : "MSS"}
    regimen_tmp = []
    if "therapeutic_regimen" in jsonDict.keys() and jsonDict["therapeutic_regimen"]:
        for regimen in jsonDict["therapeutic_regimen"]:
            if "var" in regimen.keys() and regimen["var"]:
                for a in regimen["var"]:
                    if a == mss_biomarker:
                        regimen["var"].remove(a)
            if "var" in regimen.keys() and regimen["var"]:
                regimen_tmp.append(regimen)
        jsonDict["therapeutic_regimen"] = regimen_tmp
    
    return jsonDict