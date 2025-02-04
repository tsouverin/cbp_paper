spieman = r'''
\author[a,*]{Thierry Souverin}
\author[a,b]{Jérémy Neveu}
\author[a]{Marc Betoule}
\author[a]{Sébastien Bongard}
\author[g]{Christopher W. Stubbs}
\author[g]{Elana Urbach}
\author[g]{Sasha Brownsberger}
\author[c]{Pierre Éric Blanc}
\author[e,f]{Johann Cohen-Tanugi}
\author[b]{Sylvie Dagoret-Campagne}
\author[d]{Fabrice Feinstein}
\author[a]{Delphine Hardin}
\author[a]{Claire Juramy}
\author[a]{Laurent Le Guillou}
\author[c]{Auguste Le Van Suu}
\author[b]{Marc Moniez}
\author[e]{Éric Nuss \textsuperscript{\textdagger}}
\author[e]{Bertrand Plez}
\author[a]{Nicolas Regnault}
\author[a]{Eduardo Sepulveda}
\author[e]{Kélian Sommer}
\author[x]{the LSST Dark Energy Science Collaboration}

\affil[a]{LPNHE, CNRS/IN2P3 \& Sorbonne Université, 4 place Jussieu, 75005 Paris, France}
\affil[b]{Universit\'e Paris-Saclay, CNRS, IJCLab, 91405, Orsay, France}
\affil[g]{Department of Astronomy, Harvard University, 60 Garden St., Cambridge, MA 02138, USA}
\affil[c]{Université d’Aix-Marseille \& CNRS, Observatoire de Haute-Provence, 04870 Saint Michel l’Observatoire, France}
\affil[d]{Aix Marseille Univ, CNRS/IN2P3, CPPM, Marseille, France}
\affil[e]{LUPM, Université Montpellier \& CNRS, F-34095 Montpellier, France}
\affil[f]{LPC, Université Clermont Auvergne, CNRS, F-63000 Clermont-Ferrand, France}
'''
#\affil[g]{Sorbonne Universit\'e, CNRS, Universit\'e de Paris, LPNHE, 75252 Paris Cedex 05, France}
import re

def author_format(author, institutes):
    instlist = [_i[0] for _i in institutes]
    inst = [instlist.index(_i)+1 for _i in author[0].split(',') if _i in instlist]
    return author[1] + '$^{' + ','.join([f'{_i}' for _i in inst]) + '}$'

re_institute = re.compile('affil\[(.)\]\{(.*)\}')
re_author = re.compile('author\[(.*)\]\{(.*)\}')
institutes = re_institute.findall(spieman)
authors = re_author.findall(spieman)
instlist = [_i[0] for _i in institutes]

def arxiv_format(author, institute):
    inst = [instlist.index(_i)+1 for _i in author[0].split(',') if _i in instlist]
    return author[1] + ' (' + ' and '.join([f'{_i}' for _i in inst]) + ')'
#print(institutes)
print(authors)
if True:
    with open('authors_rasti.tex', 'w') as fid:
        print('\\author{', file=fid)
        print('\n'.join([author_format(author, institutes) for author in authors]), file=fid)
        #print('}', file=fid)
        print('\n%List of institutions', file=fid)
        print(f'\n'.join([f"$^{_i+1}$ {i[1]}" for _i, i in enumerate(institutes)]), file=fid)
        print('}', file=fid)
if False:
    print(', '.join([arxiv_format(author, institutes) for author in authors]))
    print(', '.join([f'({instlist.index(_i[0])}) {_i[1]}' for _i in institutes]))
