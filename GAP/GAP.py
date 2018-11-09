#-*- coding:utf-8 -*-

import ParseTrans

def init(genomeCoorBedFile, seqFn='', showAttr=True, rem_tVersion=False, rem_gVersion=False):
    return ParseTrans.ParseTransClass(genomeCoorBedFile, seqFileName=seqFn, showAttr=showAttr, remove_tid_version=rem_tVersion, remove_gid_version=rem_gVersion)

