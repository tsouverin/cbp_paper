spieman = r'''
\author[a,*]{Thierry Souverin}
\author[b,g]{Jérémy Neveu}
\author[a]{Marc Betoule}
\author[a]{Sébastien Bongard}
\author[g]{Christopher W. Stubbs}
\author[g]{Elana Urbach}
\author[g]{Sasha Brownsberger}
\author[d]{Pierre Éric Blanc}
\author[e,f]{Johann Cohen Tanugi}
\author[b]{Sylvie Dagoret-Campagne}
\author[c]{Fabrice Feinstein}
\author[a]{Delphine Hardin}
\author[a]{Claire Juramy}
\author[a]{Laurent Le Guillou}
\author[d]{Auguste Le Van Suu}
\author[b]{Marc Moniez}
\author[e]{Bertrand Plez}
\author[a]{Nicolas Regnault}
\author[a]{Eduardo Sepulveda}
\author[e]{Kélian Sommer}



\affil[a]{LPNHE, CNRS/IN2P3 \& Sorbonne Université, 4 place Jussieu, 75005 Paris, France}
\affil[b]{Universit\'e Paris-Saclay, CNRS, IJCLab, 91405, Orsay, France}
\affil[c]{Aix Marseille Univ, CNRS/IN2P3, CPPM, Marseille, France}
\affil[d]{Université d’Aix-Marseille \& CNRS, Observatoire de Haute-Provence, 04870 Saint Michel l’Observatoire, France}
\affil[e]{LUPM, Université Montpellier \& CNRS, F-34095 Montpellier, France}
\affil[f]{LPC, université Clermont Auvergne, CNRS, F-63000 Clermont-Ferrand, France}
\affil[g]{Sorbonne Universit\'e, CNRS, Universit\'e de Paris, LPNHE, 75252 Paris Cedex 05, France}
\affil[g]]{Department of Astronomy, Harvard University, 60 Garden St., Cambridge, MA 02138, USA}
'''
import re

def author_format(author, institutes):
    instlist = [_i[0] for _i in institutes]
    inst = [instlist.index(_i)+1 for _i in author[0].split(',') if _i in instlist]
    return author[1] + '\\inst{' + ','.join([f'{_i}' for _i in inst]) + '}'

re_institute = re.compile('affil\[(.)\]\{(.*)\}')
re_author = re.compile('author\[(.*)\]\{(.*)\}')
institutes = re_institute.findall(spieman)
authors = re_author.findall(spieman)
instlist = [_i[0] for _i in institutes]

def arxiv_format(author, institute):
    inst = [instlist.index(_i)+1 for _i in author[0].split(',') if _i in instlist]
    return author[1] + ' (' + ' and '.join([f'{_i}' for _i in inst]) + ')'
#print(institutes)
#print(authors)
if True:
    print('\\author{')
    print('\n \\and '.join([author_format(author, institutes) for author in authors]))
    print('}')
    
    print('\\institute{')
    print('\n \\and '.join([_i[1] for _i in institutes]))
    print('}')
if True:
    print(', '.join([arxiv_format(author, institutes) for author in authors]))
    print(', '.join([f'({instlist.index(_i[0])}) {_i[1]}' for _i in institutes]))
