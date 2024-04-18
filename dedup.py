import xml.etree.cElementTree as ET
import hashlib as hash
import copy
import sys

def get_attr(obj,attr):
    try:
        return obj.attrib[attr]
    except KeyError:
        return None
    else:
        return None

def iter_groups(group):
    global hashlist
    global count
    rem=[]
    for child in group:
        if child.tag==rtag+'g' :#we have a group
            iter_groups(child)
        else:
            if (child.tag==rtag+'path') or (child.tag==rtag+'use'):
                h=hash.md5(str(child.attrib)).hexdigest()
                # print h
                if h in hashlist:
                    rem.append(child)
                    print "removing ",child.tag, "in",group.tag,group.attrib," -- duplicate"
                    count+=1
                else:   
                    # if (get_attr(child,'fill')!=None):
                        # if ("rgb(255,255,255)") in child.attrib['fill']:
                            # print "removing ",get_attr(child,'id'), "in",group.tag,group.attrib," -- white"
                            # rem.append(child)
                        # else:
                            # hashlist.append(h)
                    # else:
                    hashlist.append(h)
    for r in rem: 
        print "about to remove",r.attrib
        group.remove(r)
    rem=[]
    for child in group:
        if child.tag==rtag+'g' :#we have a group
            if len(child.findall('*'))==0:
                print "removing ",child.tag, "in",group.tag,group.attrib," -- empty"
                rem.append(child)
    for r in rem: group.remove(r)


def ungroup_singles(group):
    global count
    for child in group:
        #print child.tag,rtag
        if child.tag==rtag+'g' :#we have a group
            print "len(group",get_attr(child,'id'),")",len(child)
            if len(child)>1:
                ungroup_singles(child)
            else :
                if len(child)==1:
                    if (len(child[0])>=1)or(child[0].tag<>rtag+'g'):
                        print "about to promote",child[0].tag,get_attr(child[0],'id'),get_attr(child[0],'class')
                        print len(child[0])
                        moveelem=copy.deepcopy(child[0])
                        group.append(moveelem)
                        group.remove(child)
                        count+=1
                    else:
                        print "about to remove",child[0].tag,get_attr(child[0],'id'),get_attr(child[0],'class')
                    child.remove(child[0])
                    count+=1
                else:#i.e. len(child)==0
                    print "about to remove",child.tag,get_attr(child,'id'),get_attr(child,'class')
                    group.remove(child)
                    count+=1
        #else:
            # if gl==1:#and not clipped?
                #moveelem= ET.copy.deepcopy(child)

#main#
hashlist=[] 
count=0 
ET.register_namespace("","http://www.w3.org/2000/svg")
tree = ET.parse(sys.argv[1])
root = tree.getroot()
rtag= root.tag.split('}')[0]+'}'
iter_groups(root)
print "A", count," elements removed"
lcount=1
# while True:
    # count=0
    # ungroup_singles(root)
    # print lcount,":",count," empty groups removed / single elements promoted from groups"
    # lcount+=1
    # if count==0 or lcount>10:
        # break

tree.write(sys.argv[2],encoding="us-ascii", xml_declaration=True, default_namespace="", method="xml")