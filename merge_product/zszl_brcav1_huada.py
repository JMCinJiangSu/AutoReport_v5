#-*- coding:gbk -*-

def get_huada_data(raw_huada_data):
    huada_data = {}
    # 1. 检测内容：通过产品编号获取基因列表，统计基因数量
    LD0105_list = "ACVRL1, ALK, APC, ATM, ATR, ATRX, AXIN2, BAP1, BARD1, BLM, BMPR1A, BRCA1, BRCA2, BRIP1, CDC73, CDH1, CDK12, CDK4, CDKN1B, CDKN2A, CFTR, CHEK1, CHEK2, DMC1, EME1, EME2, ENG, EPAS1, EPCAM, EXT1, EXT2, FAM175A, FANCA, FANCC, FANCG, FANCI, FANCL, FH, FLCN, GALNT12, GEN1, GREM1, HOXB13, KIT, MAX, MC1R, MDH2, MEN1, MET, MITF, MLH1, MLH3, MRE11A, MSH2, MSH3, MSH6, MSR1, MUS81, MUTYH, NBN, NF1, NF2, NTHL1, NTRK1, PALB2, PDGFRA, PHOX2B, PMS1, PMS2, POLD1, POLE, PPP2R2A, PRSS1, PTCH1, PTCH2, PTEN, RAD50, RAD51B, RAD51C, RAD51D, RAD52, RAD54L, RB1, RBBP8, RET, RHBDF2, SDHA, SDHAF2, SDHB, SDHC, SDHD, SLX1A, SLX4, SMAD4, SMARCA4, SPINK1, STK11, SUFU, TMEM127, TP53, TP53BP1, TSC1, TSC2, VHL, WT1, XPC, XRCC2, XRCC3"
    LD0108_list = "ACVRL1, ALK, APC, ATM, ATRX, AXIN2, BAP1, BARD1, BLM, BMPR1A, BRCA1, BRCA2, BRIP1, CDC73, CDH1, CDK12, CDK4, CDKN1B, CDKN2A, CFTR, CHEK1, CHEK2, DMC1, EME1, EME2, ENG, EPAS1, EPCAM, EXT1, EXT2, FANCC, FANCG, FANCI, FANCL, FH, FLCN, GALNT12, GREM1, KIT, MAX, MC1R, MDH2, MEN1, MET, MITF, MLH1, MLH3, MRE11A, MSH2, MSH3, MSH6, MSR1, MUS81, MUTYH, NBN, NF1, NF2, NTHL1, NTRK1, PALB2, PDGFRA, PHOX2B, PMS1, PMS2, POLD1, POLE, PPP2R2A, PRSS1, PTCH1, PTCH2, PTEN, RAD50, RAD51B, RAD51C, RAD51D, RAD52, RAD54L, RB1, RBBP8, RET, RHBDF2, SDHA, SDHAF2, SDHB, SDHC, SDHD, SLX1A, SLX4, SMAD4, SMARCA4, SPINK1, STK11, SUFU, TMEM127, TP53, TP53BP1, TSC1, TSC2, VHL, WT1, XPC, XRCC2, XRCC3"
    gene_list = LD0105_list.split(", ") if raw_huada_data["productConfigId"] == "LD0105" else LD0108_list.split(", ") if raw_huada_data["productConfigId"] == "LD0108" else []
    huada_data["gene_count"] = len(gene_list)
    huada_data["gene_list"] = ", ".join(gene_list)
    # 2. 检测结果
    # 2.1 总结
    # 规则：读取是否报出（showFlag）为1的位点，展示3/4/5类
    # 疑问： p.?时这个字段返回什么？--> p.?时返回空
    # 疑问：拷贝数变异时怎么展示？json中用什么字段？--> a.functionAlteration in [“Duplication”, “duplication”, “DUP”, “dup”, “Deletion”, “deletion”, “DEL”, “del”] and not a.proteinResult and “c.” not in a.nucleotideResult (因为CNV和SNV的p.字段可能都为空，无法区分是SNV的重复突变还是CNV的DUP，所以再加个c.判定)
    var_sum_list = [var for var in raw_huada_data["readMutationInfo"] if var["showFlag"] == 1 and var["clinicalSignificance"] in [1, 2, 3]]
    var_sum_list = sorted(var_sum_list, key = lambda i : (i["clinicalSignificance"], i["gene"]))
    huada_data["var_sum_list"] = var_sum_list
    # 2.1 正文
    # 规则：读取是否报出（showFlag）为1、正文报出（mainTextReportFlag）为1 的3/4/5类变异
    var_main_text_list = [var for var in raw_huada_data["readMutationInfo"] if var["showFlag"] == 1 and var["mainTextReportFlag"] == 1 and var["clinicalSignificance"] in [1, 2, 3]]
    var_main_text_list = sorted(var_main_text_list, key = lambda i : (i["clinicalSignificance"], i["gene"]))
    huada_data["var_main_text_list"] = var_main_text_list
    huada_data["var_main_text_list_judge_P_LP_var"] = "yes" if [var for var in var_main_text_list if var["clinicalSignificance"] in [1, 2]] else "no"
    # 变异解析可能有重复，这边加个去重后的列表-2026.03.17
    var_main_text_inter_redup = []
    for i in var_main_text_list:
        if i["mutationAnalysis"] not in var_main_text_inter_redup:
            var_main_text_inter_redup.append(i["mutationAnalysis"])
    huada_data["var_main_text_inter_redup"] = var_main_text_inter_redup
    # 新增完成-2026.03.17
    # 3. 疾病背景简介
    # 规则：读取readT16Info，获取patientDiseaseIntroduction接口内容进行展示
    huada_data["disease_introduction"] = []
    for disease in raw_huada_data["readT16Info"]:
        disease_dict = transform_disease_introduction(disease)
        for tumor, info in disease_dict.items():
            if {"tumor": tumor, "info" : info} not in huada_data["disease_introduction"]:
                huada_data["disease_introduction"].append({"tumor": tumor, "info" : info})
    # 4.数据质控
    # 规则：读取qcResult的内容
    huada_data["qc"] = raw_huada_data["qcResult"]
    # 5. 相关基因外显子区及邻近±20bp内含子区非同义变异位点
    # 规则：读取t25和变异信息表内容的整合
    # 1. t25表读取全部内容，用“基因+c点+p点+亚区”为索引条件，在变异表中匹配，如果变异类型(变异等级)不一致，则以变异表为准
    # 2. 变异表中，比t25多的内容，且是否报出（showFlag）为1的位点，需要展示到报告中
    # 3. 根据变异等级进行排序
    # 注意！del/dup的c.，readMutationInfo中是一个空格，t25中是两个空格（最新测试数据看都只有一个空格了，还是都替换掉比较保险）
    t25_var_list = raw_huada_data["readT25Info"]
    readMutationInfo_list = [var for var in raw_huada_data["readMutationInfo"] if var["showFlag"] == 1]
    # 获取变异列表中的变异等级 
    var_clinica_dict_from_readMutationInfo = {}
    for var in readMutationInfo_list:
        key = var["gene"] + var["nucleotideResult"].replace(" ", "") + var["proteinResult"] + var["geneSubregion"]
        var_clinica_dict_from_readMutationInfo[key] = var["clinicalSignificance"]
    # 获取t25表中的变异，其中变异等级按变异列表中的为准
    clinical_dict = {
        "已知致病变异" : 1,
        "疑似致病变异" : 2,
        "意义未明变异" : 3,
        "疑似良性变异" : 4,
        "良性变异" : 5,
        "-" : 0
    }
    t25_var_key = []
    for var in t25_var_list:
        key = var["gene"] + var["nucleotideResult"].replace(" ", "") + var["proteinResult"] + var["geneSubregion"] 
        if key in var_clinica_dict_from_readMutationInfo.keys():
            var["clinicalSignificance"] = var_clinica_dict_from_readMutationInfo.get(key, "-")
        else:
            var["clinicalSignificance"] = clinical_dict.get(var["mutationClassification"], "-")
        t25_var_key.append(var["gene"] + var["nucleotideResult"].replace(" ", "") + var["proteinResult"] + var["geneSubregion"] )
        # 2026.03.17-连锁突变拆分开分别加key
        if "<br>" in var["nucleotideResult"]:
            for i in range(len(var["nucleotideResult"].split("<br>"))):
                t25_var_key.append(var["gene"] + var["nucleotideResult"].split("<br>")[i].replace(" ", "") + var["proteinResult"].split("<br>")[i] + var["geneSubregion"] )
        # 2026.03.17-完成
    permlineVariationSites_list = []
    permlineVariationSites_list.extend(t25_var_list)
    for var in readMutationInfo_list:
        key = var["gene"] + var["nucleotideResult"].replace(" ", "") + var["proteinResult"] + var["geneSubregion"] 
        if key not in t25_var_key:
            permlineVariationSites_list.append(var)
    for i in permlineVariationSites_list:
        print (i["clinicalSignificance"], i["gene"])
    huada_data["permlineVariationSites_list"] = sorted(permlineVariationSites_list, key = lambda i : (i["clinicalSignificance"], i["gene"]))

    return huada_data

# 遗传病介绍格式转化
def transform_disease_introduction(a):
    text = a["patientDiseaseIntroduction"] if "patientDiseaseIntroduction" in a.keys() and a["patientDiseaseIntroduction"] else ""
    parts = text.split("<br>")
    result = {}
    current_disease = None
    for part in parts:
        clean_part = part.strip()
        if clean_part.startswith("*") and clean_part.endswith("*") and len(clean_part) >= 2:
            current_disease = clean_part[1:-1]
            result[current_disease] = []
        else:
            if current_disease is not None:
                if part:
                    result[current_disease].append(part)             
    return result