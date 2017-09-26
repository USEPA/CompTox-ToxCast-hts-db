import os,sys
import subprocess
import shlex
import copy as cp
    
def flatten(L):
    if not L: return []
    if type(L) != list or len(L)==1: return L
    def aaa(a,b):
        if type(a)!=list: a=[a]
        if type(b)!=list: b=[b]
        return a+b

    return [i for i in list(set(reduce(aaa,L))) if i]
    
    
def ifthen(cond,if_true,if_false):
    if cond:
        return if_true
    else:
        return if_false


import smtplib
from email.mime.text import MIMEText
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart

def send_email(txt="text",subj="Job complete",email=None):
    if not email: return
    s = smtplib.SMTP('localhost')
    msg = MIMEText(txt)
    H  = dict(Subject=subj,From=email,To=email)
    for k,v in H.iteritems(): msg[k]=v
    s.sendmail(msg['From'], [msg['To']], msg.as_string())
    s.quit()

