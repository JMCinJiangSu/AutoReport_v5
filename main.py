#-*- coding:gbk -*-
import argparse
import os
import re
#import subprocess as sp
import AutoReport_Base_Code_v1_0_0

## 配置信息---------------------------------------------------------------
# 程序路径
# 代码框架程序调用方式更新，基础路径更新-2023.08.21
#BASE_DIR = "/".join((re.split("/", AutoReport_Base_Code_v1_0_0.__file__)[0:-1]))
BASE_DIR = (AutoReport_Base_Code_v1_0_0.__file__).replace("AutoReport_Base_Code_v1_0_0.py", "")

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

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

## 下面为主程序代码，用于对接报告系统，正常情况下不用动
def main(json_name, outfile):
	# 封装后linux上测试可运行，部署到报告系统会报错
	#cmd = "python %s -s %s -o %s -c %s -r %s -j %s -i %s" % (base_script, json_name, outfile, config, report_template, outjson, image)
	#sp.run(cmd, shell=True)
	# 调用方式更新-2023.08.21
	AutoReport_Base_Code_v1_0_0.get_data(json_name, outfile, config, report_template, outjson, image)

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("-s", "--json_name", dest = "json_name", required = True)
	parser.add_argument("-o", "--outfile", dest = "outfile", required = True)
	arg = parser.parse_args()

	return arg

def script_main_auto_output_report_file(json_name=None, outfile_dir=None):
	"""
	# 程序调用入口
	"""
	out_report_file, word_name = main(json_name=json_name, outfile=outfile_dir)
	return out_report_file, word_name

if __name__ == '__main__':
	args = parse_args()
	main(json_name=args.json_name, outfile=args.outfile)