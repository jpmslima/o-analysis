import re
from ete3 import Tree, TreeStyle, TextFace, NodeStyle

# 1. Carregar o arquivo Nexus
nome_arquivo = "guest_4669_seq_tree.nex"
with open(nome_arquivo, "r") as f:
    content = f.read()

# Extrair a string da árvore (usando re.DOTALL para pegar múltiplas linhas até o ';')
tree_match = re.search(r'tree TREE1\s*=\s*(.+?;)', content, re.DOTALL)
if not tree_match:
    raise ValueError("Árvore não encontrada no arquivo.")

# Limpar as quebras de linha e tabulações para o ETE3 ler como uma string contínua
newick_str = tree_match.group(1).replace('\n', '').replace('\t', '').replace('\r', '')

# Converter as anotações do TaxOnTree [&...] para o padrão NHX [&&NHX:...] 
def convert_to_nhx(match):
    inner = match.group(1)
    # HIGIENIZAÇÃO CRÍTICA: Remove aspas, parênteses (que quebram o parser) e espaços
    inner = inner.replace('"', '').replace('(', '').replace(')', '').replace(' ', '_').replace(',', ':')
    # Remove os prefixos numéricos dos nomes dos metadados (ex: 01-superkingdom -> superkingdom)
    inner = re.sub(r'(\d+)-([a-zA-Z]+)=', r'\2=', inner)
    # Remove caracteres especiais que o FigTree introduz
    inner = inner.replace('!color', 'color')
    return '[&&NHX:' + inner + ']'

newick_nhx = re.sub(r'\[\&(.*?)\]', convert_to_nhx, newick_str)

# Carregar a árvore no ETE3 (com aviso de que os nós possuem aspas no nome original)
t = Tree(newick_nhx, quoted_node_names=True)

# 2. Definir o layout para formatar os nomes das folhas
def my_layout(node):
    classe = getattr(node, "class", "")
    
    # Se for uma folha (terminal) que não foi colapsada
    if node.is_leaf():
        # Regra para Mammalia: Exibir apenas Família / Gênero
        if "Mammalia" in classe:
            # As strings estarão no formato (XX)_Nome. Pegamos a última parte.
            fam = getattr(node, "family", "Desconhecido").split("_")[-1]  
            gen = getattr(node, "genus", "Desconhecido").split("_")[-1]
            
            face = TextFace(f"{fam} / {gen}", fsize=10, fgcolor="black")
            node.add_face(face, column=0, position="branch-right")
        else:
            # Outros terminais (se houver algum solto)
            # Retiramos as aspas duplas residuais do nome original para ficar mais limpo
            face = TextFace(node.name.replace('"', ''), fsize=9)
            node.add_face(face, column=0, position="branch-right")

# 3. Aplicar as regras de colapso iterando de baixo para cima (postorder)
for node in t.traverse("postorder"):
    if node.is_leaf():
        continue
        
    superkingdom = getattr(node, "superkingdom", "")
    kingdom = getattr(node, "kingdom", "")
    classe = getattr(node, "class", "")

    # Regra A: Colapsar completamente o clado Bacteria
    if "Bacteria" in superkingdom:
        # Se o nó pai NÃO é Bacteria, significa que achamos o nó raiz de Bacteria
        parent_sk = getattr(node.up, "superkingdom", "") if node.up else ""
        if "Bacteria" not in parent_sk:
            node.img_style["draw_descendants"] = False # Executa o colapso
            node.add_face(TextFace(" Domínio Bacteria (Colapsado)", fsize=11, fgcolor="red"), column=0, position="branch-right")
            
    # Regra B: No clado Metazoa, colapsar por Classe (exceto Mammalia)
    elif "Metazoa" in kingdom and classe != "" and "Mammalia" not in classe:
        parent_class = getattr(node.up, "class", "") if node.up else ""
        # Verifica se achou o ancestral comum mais profundo daquela Classe específica
        if parent_class != classe:
            node.img_style["draw_descendants"] = False # Executa o colapso
            class_name = classe.split("_")[-1]
            # Adiciona o nome da Classe colapsada e a quantidade de folhas escondidas
            n_leaves = len(node.get_leaves())
            node.add_face(TextFace(f" Classe {class_name} ({n_leaves} spp)", fsize=11, fgcolor="blue"), column=0, position="branch-right")

# 4. Configurar estilo e renderizar a figura
ts = TreeStyle()
ts.layout_fn = my_layout
ts.show_leaf_name = False # Desativa os nomes padrão para evitar duplicação
ts.mode = "r" # "r" para retangular, "c" se quiser circular

# Gerar arquivo em PDF
t.render("arvore_filogenetica_organizada.pdf", w=200, units="mm", tree_style=ts)
print("Árvore renderizada com sucesso: arvore_filogenetica_organizada.pdf")