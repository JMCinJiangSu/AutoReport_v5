#-*- coding:gbk -*-
import os
import re
import AutoReport_Base_Code_v1_0_0
import json
import xlrd
import logging
import datetime
from bin import getTemplate

'''
    用于自动生成测试结果
    Auto_test [存放文件夹，测试集、批量生成的报告、日志文件]
        |-Auto_file.xlsx [测试集、待批量生成报告信息]
        |-log.txt [日志文件]
        |-test_report[批量生成的报告]
        |-JSON_SET [测试集json存放文件夹]
'''

Auto_dir = "./Auto_test"
if not os.path.exists(Auto_dir):
    os.mkdir(Auto_dir)

Auto_file = Auto_dir + "/Auto_file-v3.xlsx" 

test_dir = Auto_dir + "/test_report"
if not os.path.exists(test_dir):
    os.mkdir(test_dir)

logging.basicConfig(level = logging.DEBUG,
                    format = '%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s',
                    filename = Auto_dir+"/log.txt",
                    filemode = "w")

key_stran = {
	"产品名称" : "prod_name",
	"检测业务类型" : "business_type",
	"模板类型" : "report_type",
	"单位名称" : "company",
	"科室" : "hosp_depart",
	"模板名" : "report_name",
	"状态" : "status",
	"模板开发者" : "developer",
	"添加人" : "auditors",
	"添加时间" : "auditors_time",
	"备注" : "note",
	"更新记录" : "update",
	"临检" : "rummage",
	"进院" : "hospital",
	"定制" : "CustomEdition",
	"通用" : "Universal",
	"通用-简版" : "Universal_simple",
	"通用-完整" : "Universal_complete",
    "药企项目名称" : "product_name",
	"订单类型" : "order_type",
	"免费类型" : "free_type",
	"来源" : "order_origin"
	}

# 1. 获取测试数据集
# 返回格式{项目：[{样本编号:XX, json文件存放路径:XX}],}
def getTestSet():
    # 获取全部数据
    Data = []
    xls = xlrd.open_workbook(Auto_file)
    test_set_sheet = xls.sheet_by_name("test_set")
    key = test_set_sheet.row_values(0)
    for num in range(1, test_set_sheet.nrows):
        rows = test_set_sheet.row_values(num)
        tmp_dict = {}
        for i in range(len(key)):
            tmp_dict[key[i]] = rows[i]
        Data.append(tmp_dict)
    # 格式处理
    result = {}
    for i in Data:
        if int(i["statu"]) == 0:
            if i["prod_name"] not in result.keys():
                result.setdefault(i["prod_name"], [])
            result[i["prod_name"]].append(
                {
                    "sample_id" : i["sample_id"],
                    "dir_name" : i["dir_name"]
                }
            )
    return result

# 2. 获取待批量生成报告的列表  
def getTemplate_list():
    Data = []
    xls = xlrd.open_workbook(Auto_file)
    test_set_sheet = xls.sheet_by_name("template_list")
    key = test_set_sheet.row_values(0)
    for num in range(1, test_set_sheet.nrows):
        rows = test_set_sheet.row_values(num)
        tmp_dict = {}
        for i in range(len(key)):
            tmp_dict[key_stran.get(key[i])] = key_stran[rows[i]] if rows[i] in key_stran.keys() else rows[i]
        #print (tmp_dict)
        if tmp_dict["report_type"] == "Universal":
            tmp_dict["company"] = "测试医院"
        elif tmp_dict["report_type"] == "Universal_simple":
            tmp_dict["company"] = "测试医院-汇总"
        elif tmp_dict["report_type"] == "Universal_complete":
             tmp_dict["company"] = "测试医院-收费"
        tmp_dict["free_type"] = ""
        tmp_dict["order_type"] = "临检项目" if tmp_dict["order_origin"] == "LIMS" and not tmp_dict["order_type"] and tmp_dict["report_type"] == "Universal_complete" else \
                                 "汇总进院"  if tmp_dict["order_origin"] == "LIMS" and not tmp_dict["order_type"] and tmp_dict["report_type"] == "Universal_simple" else \
                                 tmp_dict["order_type"]


        if int(tmp_dict["status"]) == 0:
            Data.append(tmp_dict)

    return Data 

# 3. 创建项目-业务类型-送检单位对应测试数据
def make_test_json(temp_info):
    test_set = getTestSet().get(temp_info["prod_name"], [])
    if test_set:
        for i in test_set:
            demo_json = json.loads(open("{0}/{1}.json".format(i["dir_name"], i["sample_id"]), "r", encoding="utf-8").read())
            demo_json["sample_info"]["company"] = temp_info["company"]
            demo_json["sample_info"]["origin_company"] = temp_info["company"]
            demo_json["sample_info"]["report_module_type"] = temp_info["business_type"]
            demo_json["sample_info"]["hosp_depart"] = temp_info["hosp_depart"]
            demo_json["sample_info"]["product_name"] = temp_info["product_name"]
            demo_json["sample_info"]["free_type"] = temp_info["free_type"]
            demo_json["sample_info"]["order_type"] = temp_info["order_type"]
            demo_json["sample_info"]["order_origin"] = temp_info["order_origin"]

            # 新增药企测试
            dataJson = json.dumps(demo_json, ensure_ascii=False)
            if not os.path.exists("{0}/{1}".format(test_dir, temp_info["report_name"])):
                os.mkdir("{0}/{1}".format(test_dir, temp_info["report_name"]))
            
            outjson_dir = "{0}/{1}/{2}_{3}_{4}".format(test_dir, temp_info["report_name"], temp_info["prod_name"].replace("/", ""), \
                                                                 temp_info["business_type"], temp_info["company"]).strip()
            if not os.path.exists(outjson_dir):
                os.mkdir(outjson_dir)
            with open(outjson_dir+"/"+i["sample_id"]+".json", "w", encoding = "utf-8") as outFile:
                outFile.write(dataJson)
            i["outjson_dir"] = outjson_dir
    return test_set
    
## 配置信息---------------------------------------------------------------
# 程序路径
#BASE_DIR = "app/script/AutoReport_v5"
BASE_DIR = ""
# 配置文件
config = os.path.join(BASE_DIR, "config")
# 报告模板
report_template = os.path.join(BASE_DIR, "report_template")
# 框架代码
base_script = os.path.join(BASE_DIR, "AutoReport_Base_Code_v1_0_0.py")
# 是否要输出用于填充word的数据，T为是，其他为否
outjson = "F"
# 新增一个image文件夹
image = os.path.join(BASE_DIR, "image")
## ----------------------------------------------------------------------
# 4. 批量生成测试报告
def main():
    start = datetime.datetime.now()
    template_list = getTemplate_list()
    logging.info("**********本次共测试{0}条配置信息{1}个模板**********".format(len(template_list), len(set([i["report_name"] for i in template_list]))))
    temp_num = 1
    for i in template_list:
        logging.info("*****配置{0}/{1}测试开始*****".format(temp_num, len(template_list)))
        logging.info("      a. 项目：{0}".format(i["prod_name"]))
        logging.info("      b. 业务类型：{0}".format(i["business_type"]))
        logging.info("      c. 送检单位：{0}".format(i["company"]))
        logging.info("      d. 预计使用模板：{0}".format(i["report_name"]))
        logging.info("      e. 药企项目：{0}".format(i["product_name"] if "product_name" in i.keys() and i["product_name"] else "非药企项目"))
        logging.info("      f. 免费类型：{0}".format(i["free_type"]))
        logging.info("      g. 订单类型：{0}".format(i["order_type"] ))
        logging.info("      h. 订单来源：{0}".format(i["order_origin"]))

        test_set = make_test_json(i)
        sample_num = 1
        for sample in test_set:
            with open(sample["outjson_dir"]+"/"+sample["sample_id"]+".json", "r", encoding='utf-8') as file_json:
                jsonDict = json.load(file_json)
                report_name, merge_template, judge_brca_cnv = getTemplate.choose_template(jsonDict, config) 
                print (report_name)         
            try:
                print (sample["sample_id"], sample["outjson_dir"])
                AutoReport_Base_Code_v1_0_0.get_data(sample["sample_id"], sample["outjson_dir"], config, report_template, outjson, image)
                if report_name == i["report_name"]:
                    logging.info("样本{0}：{1}，模板调用跟预期一致，CNV类型为{2}".format(sample_num, sample["sample_id"], judge_brca_cnv))
                else:
                    logging.warning("样本{0}：{1}，模板调用跟预期不一致（{2}），注意排查！".format(sample_num, sample["sample_id"], report_name))
            except Exception as e:
                print ("ERROR")
                logging.error("样本{0}:{1}生成报告失败啦！".format(sample_num, sample["sample_id"]))
                logging.error(e)
            sample_num += 1
        logging.info("*****配置{0}/{1}测试结束*****".format(temp_num, len(template_list)))
        temp_num += 1
    end = datetime.datetime.now()
    runTime = end - start
    logging.info("**********本次测试结束啦，共耗时{0}**********".format(runTime))

if __name__ == "__main__":
    main()