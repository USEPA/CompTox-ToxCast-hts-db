import pymongo
from mongoengine import *
from util.misc import *


register_connection("hts-db","htsdb",user='devr',passwd='devr',host='localhost')

class HtsParam(EmbeddedDocument):
    meta=dict(db_alias="hts-db")
    name=StringField()
    value=DynamicField()
    desc=StringField()

class HtsAnalysisStep(Document):
    meta=dict(db_alias="hts-db")
    step=IntField()
    name=StringField(primary_key=True)
    pkg=StringField()
    fun=StringField()
    desc=StringField()
    params=ListField(EmbeddedDocumentField('HtsParam'))

class HtsAnalysisMethod(Document):
    meta=dict(db_alias="hts-db")
    name=StringField(primary_key=True)
    desc=StringField()
    steps=ListField(ReferenceField(HtsAnalysisStep))
    
class HtsChem(Document):
    meta=dict(db_alias="hts-db")
    name = StringField()
    eid  = StringField(primary_key=True)
    reps = ListField(StringField())
    chid = ObjectIdField() # txpdb Chemical
    cid  = IntField()
    casrn= StringField()
    dssid= StringField()
    ctrl = IntField()
    tags = ListField(StringField())

class HtsExp(Document):
    meta=dict(db_alias="hts-db")
    name = StringField()
    tech = StringField()
    desc = StringField()
    eid  = StringField(primary_key=True)
    src  = StringField()
    org  = StringField()
    cell = StringField()
    tags= ListField(StringField())
    ref=StringField()
    chems=ListField(ReferenceField(HtsChem))

class HtsExpLink(Document):
    meta=dict(db_alias="hts-db")
    exp  = ReferenceField(HtsExp)
    link_id= ObjectIdField()
    link_cls=StringField()
    link_db= StringField()
    tags  = ListField(StringField())
    
class HtsAssay(Document):
    meta=dict(db_alias="hts-db")
    name = StringField()
    lab  = StringField()
    desc = StringField()
    abrv = StringField()
    eid  = StringField(primary_key=True)
    txid = ListField(StringField())
    src  = StringField()
    pos_ctrl = ListField(ReferenceField('HtsChem'))
    exp  = ReferenceField(HtsExp)

class HtsResult(EmbeddedDocument):
    meta=dict(db_alias="hts-db",allow_inheritance=True)
    name=StringField()
    value=DynamicField(default=0)
    units=StringField()

class HtsConcRespNrm(EmbeddedDocument):
    meta=dict(db_alias="hts-db")
    lconc=FloatField()
    raw=FloatField()
    sm =FloatField()
    slfc=FloatField()
    ppct=FloatField()
    pzsc=FloatField()
    epct=FloatField()
    ezsc=FloatField()
    qa =FloatField()

class HtsConcRespCurveNrm(Document):
    meta=dict(db_alias="hts-db")
    chem = ReferenceField('HtsChem')
    assay = ReferenceField('HtsAssay')
    plate = ReferenceField('HtsPlate')
    wells = ListField(ReferenceField('HtsWell'))
    resp=  ListField(EmbeddedDocumentField('HtsConcRespNrm'))
    meth=  ReferenceField('HtsAnalysisMethod')
    timeh= FloatField()
    
class HtsAssayResult(Document):
    meta=dict(db_alias="hts-db")
    chem = ReferenceField('HtsChem')
    plate= ReferenceField('HtsPlate')
    assay = ReferenceField('HtsAssay')
    crc = ReferenceField('HtsConcRespCurveNrm')
    meth  = ReferenceField('HtsAnalysisMethod')
    exp  = ReferenceField('HtsExp')
    hit_call= IntField()
    lec = FloatField()
    llec = FloatField()
    lfc_eff = FloatField()
    pc_eff  = FloatField()
    llc_eff  = FloatField()
    res = ListField(EmbeddedDocumentField(HtsResult))
    

class HtsRespRaw(EmbeddedDocument):
    meta=dict(db_alias="hts-db")
    assay = ReferenceField(HtsAssay)
    raw=FloatField()
    
class HtsWell(Document):
    meta=dict(db_alias="hts-db")
    irow=IntField()
    icol=IntField()
    lconc=FloatField()
    timeh=FloatField()
    chem=ReferenceField('HtsChem')
    resp=ListField(EmbeddedDocumentField(HtsRespRaw))
    plate=ReferenceField('HtsPlate')

class HtsPlate(Document):
    meta=dict(db_alias="hts-db")
    tags = ListField(StringField)
    eid =StringField()
    timeh=FloatField()
    exp  = ReferenceField(HtsExp)
    desc=StringField()
    assays=ListField(ReferenceField(HtsAssay,dbref=True))
    meth=ReferenceField(HtsAnalysisMethod)
    rwdat = StringField()
    wells=ListField(ReferenceField(HtsWell,dbref=True))
    chems=ListField(ReferenceField(HtsChem,dbref=True))

