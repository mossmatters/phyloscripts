
import os,sys,argparse,subprocess
from ete3 import Tree,TreeStyle,TextFace,NodeStyle

helptext = '''This script will print tree figures representing the minority bipartitions
found by PhyParts. Given the species tree used for PhyParts, the PhyParts output file 'root'
and the node of interest, one PNG file will be generated per minority bipartition, and 
the decendants from the bipartition will be color-coded on the species tree.
'''

#color_it = ["Entosthodon-hungaricus-3177","Entosthodon-attenuatus-3479","Physcomitrium-spathulatum-3549","Physcomitrium-pyriforme-3727","Physcomitrium-pyriforme-3728","Physcomitrium-immersum-3176","Physcomitrium-pyriforme-3118","Physcomitrella-magdalenae-3844","Physcomitrium-hookeri-3412","Physcomitrium-sp-3842","Entosthodon-attenuatus-3835","Physcomitrium-hookeri-3409","Physcomitrium-pyriforme-3798","Physcomitrium-sp-3115","Physcomitrium-pyriforme-3387","Entosthodon-obtusus-3347","Physcomitrium-pyriforme-3404","Aphanorrhegma-serratum-3305","Physcomitridium-readeri-3892","Entosthodon-americanus-3894","Entosthodon-lindigii-3546","Physcomitrium-sp-3508","Physcomitrium-sp-3672","Entosthodon-attenuatus-3543","Physcomitrium-pyriforme-3555","Physcomitrella-patens-3403","Physcomitrium-collenchymatum-3480","Entosthodon-obtusus-3395","Physcomitrium-eurystomum-3841","Physcomitrium-sp-3551","Physcomitrium-subsphaericum-3556","Physcomitrium-collenchymatum-3178","Physcomitrium-sp-3496","Physcomitrella-patens-3139","Physcomitrium-japonicum-3413","Physcomitrium-japonicum-3411","Physcomitrium-pyriforme-3787","Physcomitrium-pyriforme-3886","Entosthodon-sp-3837","Entosthodon-subintegrus-3840","Physcomitrium-eurystomum-3392","Physcomitrium-sp-3539","Physcomitrium-pyriforme-3883","Physcomitrium-sp-3816","Entosthodon-bergianus-3509","Physcomitrium-sp-3817","Physcomitrium-sp-3814"]

def get_alternative_bipartitions(node,phyparts_root,min_alt):
    alt_biparts = []
    alt_counts = []
    for line in open(phyparts_root + ".hist.alts"):
        line = line.split()
        nodenum = int(line[2])
        if nodenum == node:
            bipart2 =  line[4].rstrip().split(",")
            bipart1 = line[3].split(":")[1].split(",")
            numalt = int(line[3].split(":")[0].replace(")","").replace("(",""))
            if numalt >= min_alt:
                alt_biparts.append((bipart1,bipart2))
                alt_counts.append(numalt)
    return alt_biparts,alt_counts
    
def render_tree(species_tree,bipart1,num_alt,png_fn,replace_taxon=None):
    color1 = "blue"
    color2 = "black"
    ts=TreeStyle()
    ts.show_leaf_name=False
    ts.show_scale=False
    nstyle = NodeStyle()
    nstyle["size"] = 0

    ts.title.add_face(TextFace("{} bipartition in {} gene trees".format(png_fn,num_alt),fsize=15,bold=True),0)
    plot_tree = species_tree
    for node in plot_tree.traverse():
        node.set_style(nstyle)
        if node.name in bipart1:
            name_face = TextFace(node.name,fgcolor=color1)
        else:
            name_face = TextFace(node.name,fgcolor=color2)
        node.add_face(name_face,0,'branch-right')
    if replace_taxon:
        for leaf in plot_tree.get_leaves:
            try:
                leaf.name=taxon_subst[leaf.name]
            except KeyError:
                continue
    plot_tree.convert_to_ultrametric()
    plot_tree.render(png_fn,tree_style=ts,w=600)        

def majority_tree(species_tree,node_num,phyparts_root):
    
    num_concord = sum([1 for line in open("{}.concord.node.{}".format(phyparts_root,node_num))])
    png_fn = "node_{}_speciestree.png".format(node_num,num_concord)
    for line in open(phyparts_root+".node.key"):
        node = int(line.split()[0])
        if node == node_num:
            subtree = Tree(line.rstrip().split()[1]+";")
            subtree_bipart = subtree.get_leaf_names()
            render_tree(species_tree,subtree_bipart,num_concord,png_fn)
                    

def main():
    parser = argparse.ArgumentParser(description=helptext,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('species_tree',help="Newick formatted species tree topology.")
    parser.add_argument('phyparts_root',help="File root name used for Phyparts.")
    parser.add_argument('node_num',type=int,default=0,help="Node number from Phyparts. To see a tree with numbered nodes, run phypartspiecharts.py with --show_nodes.")
    parser.add_argument('min_alt',type=int,default=0,help="Only print alternative bipartitions if they occur in this many gene trees")
    parser.add_argument('--taxon_subst',help="Comma-delimted file to translate tip names.")

    args = parser.parse_args()
    
    try: 
        subprocess.check_output('which convert',shell=True)
        convert = True
    except: 
        convert = False
    if args.taxon_subst:
        taxon_subst = {line.split(",")[0]:line.split(",")[1] for line in open(args.taxon_subst,'U')}
    else:
        taxon_subst = None


    alt_bipart,alt_counts = get_alternative_bipartitions(args.node_num,args.phyparts_root,args.min_alt)
    print("{} alternative bipartitions occurred in more than {} gene trees\n".format(len(alt_counts),args.min_alt))
    
    species_tree = Tree(args.species_tree)
    species_tree.ladderize(direction=1)
    majority_tree(species_tree,args.node_num,args.phyparts_root)
    
    for alt in range(len(alt_counts)):
        png_fn = "node_{}_alt_{}.png".format(args.node_num,alt)
        species_tree = Tree(args.species_tree)
        species_tree.ladderize(direction=1)
        render_tree(species_tree,alt_bipart[alt][0],alt_counts[alt],png_fn,replace_taxon=taxon_subst)
    if convert:
        os.system("convert node_{}_speciestree.png node_{}_alt_*.png node_{}.pdf".format(args.node_num,args.node_num,args.node_num))
        os.system("rm node_{}*.png".format(args.node_num))

#t = Tree("((((Discelium-nudum-3746:1,Encalypta-intermedia-3219:1)1:1.02738,((Timmia-austriaca-3619:1,Entosthodon-pulchellus-3120:1)1:0.678047,(Chamaebryum-pottioides-3630:1,Chamaebryum-pottioides-3573:1)1:5.61798)1:0.843377)1:3.51197,(((Funaria-flavicans-4092:1,Funaria-hygrometrica-3891:1)1:3.87102,((Funaria-sp-3541:1,(Funaria-hygcalvescens-3633:1,Funaria-sp-3514:1)1:0.397782)1:2.28115,(Funaria-microstoma-3834:1,(Funaria-arctica-3544:1,(Funaria-polaris-3542:1,(Funaria-arctica-3833:1,((Funaria-hygrometrica-3476:1,(Funaria-hygrometrica-3388:1,Funaria-hygrometrica-3179:1)0.38:0.00555864)1:0.739336,((Funaria-sp-3882:1,Funaria-sp-3393:1)1:0.648175,(Funaria-hygrometrica-3515:1,Funaria-hygrometrica-3632:1)0.64:0.0737125)1:0.445965)1:0.943796)1:0.575318)1:0.240042)1:0.139316)1:0.88789)1:2.00293)1:5.75785,(((Physcomitrellopsis-africana-3142:1,Entosthodon-smithhurstii-3465:1)1:2.40187,(Entosthodon-sp-3726:1,(Entosthodon-sp-3545:1,(Entosthodon-clavatus-3896:1,Entosthodon-clavatus-3895:1)1:0.669355)1:3.49935)1:0.193954)1:2.82861,(((Entosthodon-hungaricus-3177:1,Entosthodon-americanus-3894:1)1:1.55658,(Entosthodon-lindigii-3546:1,((Entosthodon-muhlenbergii-3893:1,(Entosthodon-planoconvexus-3114:1,Entosthodon-duriaei-3843:1)1:2.39883)1:3.99232,((Entosthodon-attenuatus-3835:1,(Entosthodon-attenuatus-3479:1,Entosthodon-attenuatus-3543:1)1:0.133086)1:3.23622,((Entosthodon-sp-3837:1,(Physcomitrium-sp-3842:1,Entosthodon-subintegrus-3840:1)1:2.61412)1:2.65073,(Entosthodon-bergianus-3509:1,(Entosthodon-obtusus-3395:1,Entosthodon-obtusus-3347:1)1:3.84659)0.91:0.0632156)1:0.36051)1:3.21301)0.92:0.0664078)0.81:0.0496285)1:0.110983,((Physcomitrium-hookeri-3412:1,(Physcomitrium-hookeri-3409:1,Physcomitrium-pyriforme-3404:1)1:0.375291)1:3.74685,((Physcomitridium-readeri-3892:1,(((Physcomitrium-eurystomum-3392:1,Physcomitrium-eurystomum-3841:1)0.46:0.0189278,(Physcomitrium-pyriforme-3555:1,Physcomitrium-pyriforme-3387:1)0.55:0.0260365)1:0.212244,(((Physcomitrium-sp-3816:1,Physcomitrium-sp-3508:1)1:2.44605,((Physcomitrium-sp-3551:1,Physcomitrium-sp-3539:1)1:1.2462,(Physcomitrium-japonicum-3413:1,Physcomitrium-japonicum-3411:1)0.78:0.0484692)1:1.21849)1:3.30359,((Physcomitrium-pyriforme-3118:1,Physcomitrium-pyriforme-3883:1)1:3.0447,(Physcomitrium-pyriforme-3787:1,((Physcomitrella-magdalenae-3844:1,(Physcomitrium-spathulatum-3549:1,(Physcomitrium-sp-3814:1,Physcomitrium-subsphaericum-3556:1)1:0.625116)1:0.522302)1:1.95296,(Physcomitrium-sp-3496:1,(Physcomitrium-pyriforme-3798:1,(Physcomitrium-pyriforme-3727:1,(Physcomitrium-pyriforme-3886:1,Physcomitrium-pyriforme-3728:1)1:0.129654)1:0.441593)1:1.5649)1:1.80883)0.59:0.0382)1:1.16855)1:2.09587)1:0.55998)1:0.590922)1:0.23246,(Physcomitrium-sp-3817:1,((Physcomitrella-patens-3403:1,Physcomitrella-patens-3139:1)1:4.05812,((Physcomitrium-sp-3672:1,(Aphanorrhegma-serratum-3305:1,(Physcomitrium-collenchymatum-3480:1,(Physcomitrium-sp-3115:1,Physcomitrium-collenchymatum-3178:1)0.49:0.035722)1:2.5141)1:0.389341)0.88:0.152428,Physcomitrium-immersum-3176:1)1:0.479012)1:1.1982)1:2.30257)1:0.817867)1:1.25824)1:4.29516)1:0.799009)1:5.4761)1:1,Goniomitrium-africanum-4081:1);")
#t.ladderize(direction=1)
#ts = TreeStyle()
#ts.show_leaf_name = False

if __name__ == "__main__":main()
