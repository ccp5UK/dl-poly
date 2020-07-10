import sys
import re
from collections import deque
import logging

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

R_USE = 0
R_VAR = 0
VAR_RE = re.compile(r" *(?P<var>[a-zA-Z_0-9]+) *(?P<rest>(?:\((?P<param>(?:[^()]+|\((?:[^()]+|\([^()]*\))*\))*)\))? *(?:= *(?P<value>(:?[^\"',()\[\]]+|\[(?:[^()\"']+|\([^()\"']*\)|\"[^\"]*\"|'[^']*')*\]|\((?:[^()\"']+|\([^()\"']*\)|\"[^\"]*\"|'[^']*')*\)|\"[^\"]*\"|'[^']*')+))?)? *(?:(?P<continue>,)|\n?) *", re.IGNORECASE)
USE_PARSE_RE = re.compile(
        r" *use( +|(?P<intrinsic> *, *Intrinsic * :: *))(?P<module>[a-zA-Z_][a-zA-Z_0-9]*)(?P<only> *, *only *:)? *(?P<imports>.*)$",
    flags=re.IGNORECASE)
commonUsesRe = re.compile(
    "^#include *\"([^\"]*(cp_common_uses.f90|base_uses.f90))\"")
localNameRe = re.compile(
    " *(?P<localName>[a-zA-Z_0-9]+)(?: *= *> *[a-zA-Z_0-9]+)? *$")
operatormRe = re.compile("\s*(?P<localName>Operator\s*\(.*\))?\s*$",flags=re.I)
VAR_DECL_RE = re.compile(
    r" *(?P<type>integer(?: *\* *[0-9]+)?|logical|character(?: *\* *[0-9]+)?|real(?: *\* *[0-9]+)?|complex(?: *\**[0-9]+)?|type|class) *(?P<parameters>\((?:[^()]+|\((?:[^()]+|\([^()]*\))*\))*\))? *(?P<attributes>(?: *, *[a-zA-Z_0-9]+(?: *\((?:[^()]+|\((?:[^()]+|\([^()]*\))*\))*\))?)+)? *(?P<dpnt>::)?(?P<vars>[^\n]+)\n?", re.IGNORECASE)  # $
INDENT_SIZE = 2
DECL_LINELENGTH = 100

OMP_DIR_RE = re.compile(r"^\s*(!\$omp)", re.IGNORECASE)
OMP_RE = re.compile(r"^\s*(!\$)", re.IGNORECASE)


class CharFilter(object):
    """
    An iterator to wrap the iterator returned by `enumerate`
    and ignore comments and characters inside strings
    """

    def __init__(self, it):
        self._it = it
        self._instring = ''

    def __iter__(self):
        return self

    def __next__(self):
        """ python 3 version"""
        pos, char = next(self._it)
        if not self._instring and char == '!':
            raise StopIteration

        # detect start/end of a string
        if char == '"' or char == "'":
            if self._instring == char:
                self._instring = ''
            elif not self._instring:
                self._instring = char

        if self._instring:
            return self.__next__()

        return (pos, char)

    def next(self):
        """ python 2 version"""
        pos, char = self._it.next()
        if not self._instring and char == '!':
            raise StopIteration

        # detect start/end of a string
        if char == '"' or char == "'":
            if self._instring == char:
                self._instring = ''
            elif not self._instring:
                self._instring = char

        if self._instring:
            return self.next()

        return (pos, char)


class InputStream(object):
    """
    Class to read logical Fortran lines from a Fortran file.
    """

    def __init__(self, infile):
        self.line_buffer = deque([])
        self.infile = infile
        self.line_nr = 0

    def nextFortranLine(self):
        """Reads a group of connected lines (connected with &, separated by newline or semicolon)
        returns a touple with the joined line, and a list with the original lines.
        Doesn't support multiline character constants!
        """
        lineRe = re.compile(
            # $
            r"(?:(?P<preprocessor>#.*\n?)| *(&)?(?P<core>(?:!\$|[^&!\"']+|\"[^\"]*\"|'[^']*')*)(?P<continue>&)? *(?P<comment>!.*)?\n?)",
            re.IGNORECASE)
        joinedLine = ""
        comments = []
        lines = []
        continuation = 0

        while 1:
            if not self.line_buffer:
                line = self.infile.readline().replace("\t", 8 * " ")
                self.line_nr += 1
                # convert OMP-conditional fortran statements into normal fortran statements
                # but remember to convert them back
                is_omp_conditional = False
                omp_indent = 0
                if OMP_RE.match(line):
                    omp_indent = len(line) - len(line.lstrip(' '))
                    line = OMP_RE.sub('', line, count=1)
                    is_omp_conditional = True
                line_start = 0
                for pos, char in CharFilter(enumerate(line)):
                    if char == ';' or pos + 1 == len(line):
                        self.line_buffer.append(omp_indent * ' ' + '!$' * is_omp_conditional +
                                                line[line_start:pos + 1])
                        omp_indent = 0
                        is_omp_conditional = False
                        line_start = pos + 1
                if(line_start < len(line)):
                    # line + comment
                    self.line_buffer.append('!$' * is_omp_conditional +
                                            line[line_start:])

            if self.line_buffer:
                line = self.line_buffer.popleft()

            if not line:
                break

            lines.append(line)
            m = lineRe.match(line)
            if not m or m.span()[1] != len(line):
                # FIXME: does not handle line continuation of
                # omp conditional fortran statements
                # starting with an ampersand.
                raise SyntaxError("unexpected line format:" + repr(line))
            if m.group("preprocessor"):
                if len(lines) > 1:
                    raise SyntaxError(
                        "continuation to a preprocessor line not supported " + repr(line))
                comments.append(line)
                break
            coreAtt = m.group("core")
            if OMP_RE.match(coreAtt) and joinedLine.strip():
                # remove omp '!$' for line continuation
                coreAtt = OMP_RE.sub('', coreAtt, count=1).lstrip()
            joinedLine = joinedLine.rstrip("\n") + coreAtt
            if coreAtt and not coreAtt.isspace():
                continuation = 0
            if m.group("continue"):
                continuation = 1
            if line.lstrip().startswith('!') and not OMP_RE.search(line):
                comments.append(line.rstrip('\n'))
            elif m.group("comment"):
                comments.append(m.group("comment"))
            else:
                comments.append('')
            if not continuation:
                break
        return (joinedLine, comments, lines)


def parseRoutine(inFile):
    """Parses a routine"""
    logger = logging.getLogger('prettify-logger')

    FCT_RE = re.compile(
        r"^([^\"'!]* )?FUNCTION\s+\w+\s*(\(.*\))?(\s*RESULT\s*\(\w+\))?\s*;?\s*$",
        re.IGNORECASE)

    SUBR_RE = re.compile(
        r"^([^\"'!]* )?SUBROUTINE\s+\w+\s*(\(.*\))?\s*;?\s*$", re.IGNORECASE)

    endRe = re.compile(r" *end\s*(?:subroutine|function)", re.IGNORECASE)
    startRoutineRe = re.compile(r"^([^\"'!]* )?(?P<kind>subroutine|function) +(?P<name>[a-zA-Z_][a-zA-Z_0-9]*) *(?:\((?P<arguments>[^()]*)\))? *(?:result *\( *(?P<result>[a-zA-Z_][a-zA-Z_0-9]*) *\))? *(?:bind *\([^()]+\))? *\n?", re.IGNORECASE)  # $
    typeBeginRe = re.compile(r" *(?P<type>integer(?: *\* *[0-9]+)?|logical|character(?: *\* *[0-9]+)?|real(?: *\**[0-9]+)?|complex(?: *\* *[0-9]+)?|type|class)[,( ]",
                             re.IGNORECASE)
    attributeRe = re.compile(
        r" *, *(?P<attribute>[a-zA-Z_0-9]+) *(?:\( *(?P<param>(?:[^()]+|\((?:[^()]+|\([^()]*\))*\))*)\))? *", re.IGNORECASE)
    ignoreRe = re.compile(r" *(?:|implicit +none *)$", re.IGNORECASE)
    interfaceStartRe = re.compile(r" *interface *$", re.IGNORECASE)
    interfaceEndRe = re.compile(r" *end +interface *$", re.IGNORECASE)
    routine = {'preRoutine': [],
               'core': [],
               'strippedCore': [],
               'begin': [],
               'end': [],
               'preDeclComments': [],
               'declarations': [],
               'declComments': [],
               'postDeclComments': [],
               'parsedDeclarations': [],
               'postRoutine': [],
               'kind': None, 'name': None, 'arguments': None, 'result': None,
               'interfaceCount': 0,
               'use': []
               }
    includeRe = re.compile(
        r"#? *include +[\"'](?P<file>.+)[\"'] *$", re.IGNORECASE)
    stream = InputStream(inFile)
    while 1:
        (jline, _, lines) = stream.nextFortranLine()
        if len(lines) == 0:
            break
        if FCT_RE.match(jline) or SUBR_RE.match(jline):
            break
        routine['preRoutine'].extend(lines)
        m = includeRe.match(lines[0])
        if m:
            try:
                subF = open(m.group('file'), 'r')
                subStream = InputStream(subF)
                while 1:
                    (subjline, _, sublines) = subStream.nextFortranLine()
                    if not sublines:
                        break
                    routine['strippedCore'].append(subjline)
                subF.close()
            except:
                import traceback
                logger.debug(
                    "error trying to follow include " + m.group('file') + '\n')
                logger.debug(
                    "warning this might lead to the removal of used variables\n")
                if logger.isEnabledFor(logging.DEBUG):
                    traceback.print_exc()
    if jline:
        routine['begin'] = lines
        m = startRoutineRe.match(jline)
        if not m or m.span()[1] != len(jline):
            raise SyntaxError(
                "unexpected subroutine start format:" + repr(lines))
        routine['name'] = m.group('name')
        routine['kind'] = m.group('kind')
        if (m.group('arguments') and m.group('arguments').strip()):
            routine['arguments'] = list(map(lambda x: x.strip(),
                                            m.group('arguments').split(",")))
        if (m.group('result')):
            routine['result'] = m.group('result')
        if (not routine['result'])and(routine['kind'].lower() == "function"):
            routine['result'] = routine['name']
    while 1:
        (jline, comment_list, lines) = stream.nextFortranLine()
        comments = '\n'.join(_ for _ in comment_list)
        if len(lines) == 0:
            break
        if lines[0].lower().startswith("#include"):
            break
        if not ignoreRe.match(jline):
            if typeBeginRe.match(jline):
                if routine['postDeclComments']:
                    routine['declComments'].extend(routine['postDeclComments'])
                routine['postDeclComments'] = []

            if typeBeginRe.match(jline):
                m = VAR_DECL_RE.match(jline)
                if (m.group('type').lower() == 'type' and
                        not m.group('parameters')):
                    break
                if not m or m.span()[1] != len(jline):
                    raise SyntaxError("unexpected type format:" + repr(jline))
                decl = {'type': m.group("type"),
                        'parameters': None,
                        'attributes': [],
                        'vars': []}
                if m.group('parameters'):
                    decl['parameters'] = (m.group("parameters").replace(" ", "").
                                          replace(",", ", "))
                str = m.group("attributes")
                while(str):
                    m2 = attributeRe.match(str)
                    if not m2:
                        raise SyntaxError("unexpected attribute format " +
                                          repr(str) + " in " + repr(lines))
                    decl['attributes'].append(m2.group().replace(" ", "").
                                              replace(",", ", ")[2:])
                    str = str[m2.span()[1]:]
                str = m.group("vars")
                while 1:
                    m2 = VAR_RE.match(str)
                    if not m2:
                        raise SyntaxError("unexpected var format " +
                                          repr(str) + " in " + repr(lines))
                    var = m2.group("var")
                    if m2.group("param"):
                        var += "(" + m2.group("param") + ")"
                    if m2.group("value"):
                        var += " = "
                        var += m2.group("value")
                    decl['vars'].append(var)
                    str = str[m2.span()[1]:]
                    if not m2.group("continue"):
                        if str:
                            raise SyntaxError("error parsing vars (leftover=" +
                                              repr(str) + ") in " + repr(lines))
                        break
                routine['parsedDeclarations'].append(decl)
            elif interfaceStartRe.match(jline):
                istart = lines
                interfaceDeclFile = StringIO()
                while 1:
                    (jline, _, lines) = stream.nextFortranLine()
                    if interfaceEndRe.match(jline):
                        iend = lines
                        break
                    interfaceDeclFile.writelines(lines)
                interfaceDeclFile = StringIO(interfaceDeclFile.getvalue())
                iroutines = []
                while 1:
                    iroutine = parseRoutine(interfaceDeclFile)
                    if not iroutine['kind']:
                        if len(iroutines) == 0:
                            interfaceDeclFile.seek(0)
                            raise SyntaxError("error parsing interface:" +
                                              repr(interfaceDeclFile.read()))
                        iroutines[-1]['postRoutine'].extend(
                            iroutine['preRoutine'])
                        break
                    iroutines.append(iroutine)
                for iroutine in iroutines:
                    routine['interfaceCount'] += 1
                    decl = {'type': 'z_interface%02d' % (routine['interfaceCount']),
                            'parameters': None,
                            'attributes': [],
                            'vars': [iroutine['name']],
                            'iroutine': iroutine,
                            'istart': istart,
                            'iend': iend
                            }
                    routine['parsedDeclarations'].append(decl)
            elif USE_PARSE_RE.match(jline):
                routine['use'].append("".join(lines))
            else:
                break
        routine['declarations'].append("".join(lines))
        if (len(routine['parsedDeclarations']) == 0 and len(routine['use']) == 0 and
                not re.match(" *implicit +none *$", jline, re.IGNORECASE)):
            routine['preDeclComments'].append("".join(lines))
        else:
            routine['postDeclComments'].append(comments)
    containsRe = re.compile(r" *contains *$", re.IGNORECASE)

    while len(lines) > 0:
        if endRe.match(jline):
            routine['end'] = lines
            break
        routine['strippedCore'].append(jline)
        routine['core'].append("".join(lines))
        if containsRe.match(lines[0]):
            break
        m = includeRe.match(lines[0])
        if m:
            try:
                subF = open(m.group('file'), 'r')
                subStream = InputStream(subF)
                while 1:
                    (subjline, _, sublines) = subStream.nextFortranLine()
                    if not sublines:
                        break
                    routine['strippedCore'].append(subjline)
                subF.close()
            except:
                import traceback
                logger.debug(
                    "error trying to follow include " + m.group('file') + '\n')
                logger.debug(
                    "warning this might lead to the removal of used variables\n")
                if logger.isEnabledFor(logging.DEBUG):
                    traceback.print_exc()
        (jline, _, lines) = stream.nextFortranLine()
    return routine


def findWord(word, text, options=re.IGNORECASE):
    """Returns the position of word in text or -1 if not found.
    A match is valid only if it is a whole word (i.e. findWord('try','retry')
    returns false)"""
    wordRe = re.compile("(?<![a-zA-Z_0-9%])" + word +
                        "(?![a-zA-Z_0-9])|(?<=[0-9.]_)" + word + "(?![a-zA-Z_0-9])", options)
    m = wordRe.search(text)
    if m:
        pos = m.span()[0]
    else:
        pos = -1
    return pos


def enforceDeclDependecies(declarations):
    """enforces the dependencies between the vars
    and compacts the declarations, returns the variables needed by other variables"""
    idecl = 0
    ii = 0
    while idecl < len(declarations):
        typeParam = "".join(declarations[idecl]['attributes'])
        if declarations[idecl]['parameters']:
            typeParam += " " + declarations[idecl]['parameters']
        typeParam = typeParam.lower()

        ivar = 0
        while ivar < len(declarations[idecl]['vars']):
            moved = 0
            m = VAR_RE.match(declarations[idecl]['vars'][ivar])
            if not m:
                raise SyntaxError('could not match var ' +
                                  repr(declarations[idecl]['vars'][ivar]))
            rest = m.group("rest")
            rest = rest.lower()
            if rest:
                for ivar2 in range(ivar + 1, len(declarations[idecl]['vars'])):
                    m = VAR_RE.match(declarations[idecl]['vars'][ivar2])
                    if findWord(m.group('var').lower(), rest) != -1:
                        moved = ivar2 + 1
            if moved:
                declarations[idecl]['vars'][moved:moved] = [
                    declarations[idecl]['vars'][ivar]]
                del declarations[idecl]['vars'][ivar]
            else:
                for idecl2 in range(idecl + 1, len(declarations)):
                    for ivar2 in range(len(declarations[idecl2]['vars'])):
                        ii += 1
                        if ii > 100000:
                            raise Error("could not enforce all constraints")
                        m = VAR_RE.match(declarations[idecl2]['vars'][ivar2])
                        if (ivar == 0 and
                                findWord(m.group('var').lower(), typeParam) != -1):
                            declarations.insert(
                                idecl2 + 1, declarations[idecl])
                            del declarations[idecl]
                            ivar = 0
                            moved = 1
                            break
                        if rest and findWord(m.group('var').lower(), rest) != -1:
                            if len(declarations[idecl]['vars']) > 1:
                                newDecl = {}
                                newDecl.update(declarations[idecl])
                                newDecl['vars'] = [
                                    declarations[idecl]['vars'][ivar]]
                                declarations.insert(idecl2 + 1, newDecl)
                                del declarations[idecl]['vars'][ivar]
                            else:
                                declarations.insert(idecl2 + 1,
                                                    declarations[idecl])
                                del declarations[idecl]
                                ivar = 0
                            moved = 1
                            break
                    if moved:
                        break
            if not moved:
                ivar += 1
        idecl += 1

    for i in range(len(declarations) - 1, 0, -1):
        if (declarations[i]['normalizedType'].lower() ==
                declarations[i - 1]['normalizedType'].lower()):
            declarations[i - 1]['vars'].extend(declarations[i]['vars'])
            del declarations[i]


def sortDeclarations(declarations):
    """sorts, compacts declarations and respects dependencies
    normalizedType has to be defined for the declarations"""

    declarations.sort(key=lambda x: x['normalizedType'].lower())

    for i in range(len(declarations) - 1, 0, -1):
        if (declarations[i]['normalizedType'].lower() ==
                declarations[i - 1]['normalizedType'].lower()):
            declarations[i - 1]['vars'].extend(declarations[i]['vars'])
            del declarations[i]

    for decl in declarations:
        decl['vars'].sort(key=lambda x: x.lower())
    enforceDeclDependecies(declarations)


def writeRoutine(routine, outFile):
    """writes the given routine to outFile"""
    outFile.writelines(routine["preRoutine"])
    outFile.writelines(routine["begin"])
    outFile.writelines(routine["declarations"])
    outFile.writelines(routine["core"])
    outFile.writelines(routine["end"])
    outFile.writelines(routine["postRoutine"])


def writeInCols(dLine, indentCol, maxCol, indentAtt, file):
    """writes out the strings (trying not to cut them) in dLine up to maxCol
    indenting each newline with indentCol.
    The '&' of the continuation line is at maxCol.
    indentAtt is the actual intent, and the new indent is returned"""

    strRe = re.compile(r"('[^'\n]*'|\"[^\"\n]*\")")
    nonWordRe = re.compile(r"(\(/|/\)|[^-+a-zA-Z0-9_.])")
    maxSize = maxCol - indentCol - 1
    tol = min(maxSize / 6, 6) + indentCol
    for fragment in dLine:
        if indentAtt + len(fragment) < maxCol:
            file.write(fragment)
            indentAtt += len(fragment)
        elif len(fragment.lstrip()) <= maxSize:
            file.write("&\n" + (" " * indentCol))
            file.write(fragment.lstrip())
            indentAtt = indentCol + len(fragment.lstrip())
        else:
            sPieces = strRe.split(fragment)
            for sPiece in sPieces:
                if sPiece and (not (sPiece[0] == '"' or sPiece[0] == "'")):
                    subPieces = nonWordRe.split(sPiece)
                else:
                    subPieces = [sPiece]
                for subPiece in subPieces:
                    if indentAtt == indentCol:
                        file.write(subPiece.lstrip())
                        indentAtt += len(subPiece.lstrip())
                    elif indentAtt < tol or indentAtt + len(subPiece) < maxCol:
                        file.write(subPiece)
                        indentAtt += len(subPiece)
                    else:
                        file.write("&\n" + (" " * indentCol))
                        file.write(subPiece.lstrip())
                        indentAtt = indentCol + len(subPiece.lstrip())
    return indentAtt


def writeCompactDeclaration(declaration, file,maxoffset):
    """Writes a declaration in a compact way"""
    d = declaration
    if 'iroutine' in d.keys():
        file.writelines(d['istart'])
        writeRoutine(d['iroutine'], file)
        file.writelines(d['iend'])
    else:
        if len(d['vars']) > 0:
            decl = " " * INDENT_SIZE * 2 + d['type']
            if d['parameters']:  # do not drop empty parameter lists?
                decl += d['parameters']
            if d['attributes']:
                for a in d['attributes']:
                    decl += ", " + a
            decl += " :: "

            dLine = [decl]
            for var in d['vars']:
                cur_len = sum([len(l) for l in dLine])
                if(len(dLine) > 1 and cur_len + len(var) > 600):
                    writeInCols(dLine, 3 * INDENT_SIZE,
                                DECL_LINELENGTH, 0, file)
                    file.write("\n")
                    dLine = [decl]
                if(len(dLine) > 1):
                    dLine[-1] += ", "
                dLine.append(var)
            writeInCols(dLine, 3 * INDENT_SIZE, DECL_LINELENGTH, 0, file)
            file.write("\n")


def writeExtendedDeclaration(declaration, file,maxoffset):
    """Writes a declaration in a nicer way (using more space)"""
    d = declaration
    if len(d['vars']) == 0:
        return
    if 'iroutine' in d.keys():
        file.writelines(d['istart'])
        writeRoutine(d['iroutine'], file)
        file.writelines(d['iend'])
    else:
        dLine = []
        dLine.append(" " * INDENT_SIZE * 2 + d['type'])
        if d['parameters']:  # do not drop empty parameter lists?
            dLine.append(d['parameters'])
        if d['attributes']:
            for a in d['attributes']:
                dLine[-1:] = [dLine[-1] + ", "]
                dLine.append(a)
        indentAtt = writeInCols(dLine, 3 * INDENT_SIZE,
                                maxoffset + 1 + 2 * INDENT_SIZE, 0, file)
        file.write(" " * (maxoffset + 2 * INDENT_SIZE - indentAtt))
        file.write(" :: ")
        indentAtt = maxoffset + 8

        dLine = []
        for var in d['vars'][:-1]:
            dLine.append(var + ", ")
        dLine.append(d['vars'][-1])
        writeInCols(dLine, maxoffset + 4 + 2 * INDENT_SIZE,
                    DECL_LINELENGTH, indentAtt, file)
        file.write("\n")

def writeExtendedDeclarationA(declaration, file,offset):
    """Writes a declaration in a nicer way (using more space)"""
    d = declaration
    if len(d['vars']) == 0:
        return
    if 'iroutine' in d.keys():
        file.writelines(d['istart'])
        writeRoutine(d['iroutine'], file)
        file.writelines(d['iend'])
    else:
        dLine = []
        dLine.append(" " * INDENT_SIZE * 2 + d['type'])
        if d['parameters']:  # do not drop empty parameter lists?
            dLine.append(d['parameters'])
        intent = False
        i = 0
        if d['attributes']:
            sintent = ''
            for a in d['attributes']:
                if a[0:6] != 'Intent':
                   dLine[-1:] = [dLine[-1] + ", "]
                   dLine.append(a)
                else:
                   intent = True
                   sintent = a

        indentAtt = writeInCols(dLine, 3 * INDENT_SIZE,
                                offset + 1 + 2 * INDENT_SIZE, 0, file)
        if intent:
            file.write(", ")
            i = 15
        file.write(" " * (offset + 2 * INDENT_SIZE - indentAtt-i))
        if intent:
            a = re.search('\(([^\)]+)\)',sintent,re.I)
            m = a.group(1)
            if (len(m) ==  2 ):
                sintent = "Intent(In   )"
            if (len(m) ==  3 ):
                sintent = "Intent(  Out)"
            file.write(sintent)
        file.write(" :: ")
        indentAtt = offset + 8

        dLine = []
        for var in d['vars'][:-1]:
            dLine.append(var + ", ")
        dLine.append(d['vars'][-1])
        writeInCols(dLine, offset + 4 + 2 * INDENT_SIZE,
                    DECL_LINELENGTH, indentAtt, file)
        file.write("\n")

def writeDeclarationsA(parsedDeclarations, file):
    """Writes the declarations to the given file"""
    maxoffset = 0
    for d in parsedDeclarations:
        l = len(d['normalizedType'])
        a = re.search('Intent\(([^\)]+)\)',d['normalizedType'],re.I)
        if a:
            m = a.group(1)
            if len(m) == 2:
                l = l + 3
            if len(m) == 3:
                l = l + 2
        maxoffset = max(maxoffset,l)
    for d in parsedDeclarations:
        maxLenVar = 0
        totalLen = 0
        for v in d['vars']:
            maxLenVar = max(maxLenVar, len(v))
            totalLen += len(v)
        #if maxLenVar > 30 or totalLen > DECL_LINELENGTH - 4:
        #    writeCompactDeclarationA(d, file)
        #else:
        #    writeExtendedDeclarationA(d, file,maxoffset)
        writeExtendedDeclarationA(d, file,maxoffset)



def writeDeclarations(parsedDeclarations, file):
    """Writes the declarations to the given file"""
    maxoffset = 0
    for d in parsedDeclarations:
        maxoffset = max(maxoffset,len(d['normalizedType']))
    for d in parsedDeclarations:
        maxLenVar = 0
        totalLen = 0
        for v in d['vars']:
            maxLenVar = max(maxLenVar, len(v))
            totalLen += len(v)
        #if maxLenVar > 30 or totalLen > DECL_LINELENGTH - 4:
        #    writeCompactDeclaration(d, file,maxoffset)
        #else:
        #    writeExtendedDeclaration(d, file,maxoffset)
        writeExtendedDeclaration(d, file,maxoffset)


def cleanDeclarations(routine):
    """cleans up the declaration part of the given parsed routine
    removes unused variables"""
    logger = logging.getLogger('prettify-logger')

    global R_VAR
    containsRe = re.compile(r" *contains *$", re.IGNORECASE)
    if routine['core']:
        if containsRe.match(routine['core'][-1]):
            logger.debug("routine %s contains other routines\ndeclarations not cleaned\n" %
                         (routine['name']))
            return
    commentToRemoveRe = re.compile(
        r" *! *(?:interface|arguments|parameters|locals?|\** *local +variables *\**|\** *local +parameters *\**) *$", re.IGNORECASE)
    nullifyRe = re.compile(
        r" *nullify *\(([^()]+)\) *\n?", re.IGNORECASE | re.MULTILINE)

    if not routine['kind']:
        return
    if (routine['core']):
        if re.match(" *type *[a-zA-Z_]+ *$", routine['core'][0], re.IGNORECASE):
            logger.debug("routine %s contains local types, not fully cleaned\n" %
                         (routine['name']))
        if re.match(" *import+ *$", routine['core'][0], re.IGNORECASE):
            logger.debug("routine %s contains import, not fully cleaned\n" %
                         (routine['name']))
    if re.search("^#", "".join(routine['declarations']), re.MULTILINE):
        logger.debug("routine %s declarations contain preprocessor directives\ndeclarations not cleaned\n" % (
            routine['name']))
        return
    try:
        rest = "".join(routine['strippedCore']).lower()
        nullifys = ",".join(nullifyRe.findall(rest))
        rest = nullifyRe.sub("", rest)
        paramDecl = []
        decls = []
        for d in routine['parsedDeclarations']:
            d['normalizedType'] = d['type']
            if d['parameters']:
                d['normalizedType'] += d['parameters']
            if (d["attributes"]):
                d['attributes'].sort(key=lambda x: x.lower())
                d['normalizedType'] += ', '
                d['normalizedType'] += ', '.join(d['attributes'])
            if "parameter" in map(str.lower, d['attributes']):
                paramDecl.append(d)
            else:
                decls.append(d)

        sortDeclarations(paramDecl)
        sortDeclarations(decls)
        has_routinen = 0
        pos_routinep = -1
        for d in paramDecl:
            for i in range(len(d['vars'])):
                v = d['vars'][i]
                m = VAR_RE.match(v)
                lowerV = m.group("var").lower()
                if lowerV == "routinen":
                    has_routinen = 1
                    d['vars'][i] = "routineN = '" + routine['name'] + "'"
                elif lowerV == "routinep":
                    pos_routinep = i
                    d['vars'][i] = "routineP = moduleN//':'//routineN"
            if not has_routinen and pos_routinep >= 0:
                d['vars'].insert(
                    pos_routinep, "routineN = '" + routine['name'] + "'")

        if routine['arguments']:
            routine['lowercaseArguments'] = list(map(
                lambda x: x.lower(), routine['arguments']))
        else:
            routine['lowercaseArguments'] = []
        if routine['result']:
            routine['lowercaseArguments'].append(routine['result'].lower())
        argDeclDict = {}
        localDecl = []
        for d in decls:
            localD = {}
            localD.update(d)
            localD['vars'] = []
            argD = None
            for v in d['vars']:
                m = VAR_RE.match(v)
                lowerV = m.group("var").lower()
                if lowerV in routine['lowercaseArguments']:
                    argD = {}
                    argD.update(d)
                    argD['vars'] = [v]
                    if lowerV in argDeclDict.keys():
                        raise SyntaxError(
                            "multiple declarations not supported. var=" + v +
                            " declaration=" + str(d) + "routine=" + routine['name'])
                    argDeclDict[lowerV] = argD
                else:
                    pos = findWord(lowerV, rest)
                    if (pos != -1):
                        localD['vars'].append(v)
                    else:
                        if findWord(lowerV, nullifys) != -1:
                            if not rmNullify(lowerV, routine['core']):
                                raise SyntaxError(
                                    "could not remove nullify of " + lowerV +
                                    " as expected, routine=" + routine['name'])
                        logger.info("removed var %s in routine %s\n" %
                                    (lowerV, routine['name']))
                        R_VAR += 1
            if (len(localD['vars'])):
                localDecl.append(localD)
        argDecl = []
        for arg in routine['lowercaseArguments']:
            if arg in argDeclDict.keys():
                argDecl.append(argDeclDict[arg])
            else:
                logger.debug("warning, implicitly typed argument '" +
                             arg + "' in routine " + routine['name'] + '\n')
        if routine['kind'].lower() == 'function':
            aDecl = argDecl[:-1]
        else:
            aDecl = argDecl

        # try to have arg/param/local, but checks for dependencies arg/param
        # and param/local
        argDecl.extend(paramDecl)
        enforceDeclDependecies(argDecl)
        splitPos = 0
        for i in range(len(argDecl) - 1, -1, -1):
            if not 'parameter' in map(str.lower, argDecl[i]['attributes']):
                splitPos = i + 1
                break
        paramDecl = argDecl[splitPos:]
        argDecl = argDecl[:splitPos]
        paramDecl.extend(localDecl)
        enforceDeclDependecies(paramDecl)
        splitPos = 0
        for i in range(len(paramDecl) - 1, -1, -1):
            if 'parameter' in map(str.lower, paramDecl[i]['attributes']):
                splitPos = i + 1
                break
        localDecl = paramDecl[splitPos:]
        paramDecl = paramDecl[:splitPos]

        newDecl = StringIO()
        for comment in routine['preDeclComments']:
            if not commentToRemoveRe.match(comment):
                newDecl.write(comment)
        newDecl.writelines(routine['use'])
        writeDeclarationsA(argDecl, newDecl)
        if argDecl and paramDecl:
            newDecl.write("\n")
        writeDeclarations(paramDecl, newDecl)
        if (argDecl or paramDecl) and localDecl:
            newDecl.write("\n")
        writeDeclarations(localDecl, newDecl)
        if argDecl or paramDecl or localDecl:
            newDecl.write("\n")
        wrote = 0
        for comment in routine['declComments']:
            if comment.strip() and not commentToRemoveRe.match(comment):
                newDecl.write(comment.strip())
                newDecl.write("\n")
                wrote = 1
        if wrote:
            newDecl.write("\n")
        routine['declarations'] = [newDecl.getvalue()]
    except:
        if 'name' in routine.keys():
            logger.critical("exception cleaning routine " +
                            routine['name'])
        logger.critical("parsedDeclartions=" +
                        str(routine['parsedDeclarations']))
        raise

    newDecl = StringIO()
    if routine['postDeclComments']:
        comment_start = 0
        for comment in routine['postDeclComments']:
            if comment.strip():
                break
            else:
                comment_start += 1

        for comment in routine['postDeclComments'][comment_start:]:
            if not commentToRemoveRe.match(comment):
                newDecl.write(comment)
                newDecl.write("\n")
        routine['declarations'][0] += newDecl.getvalue()


def rmNullify(var, strings):
    removed = 0
    var = var.lower()
    nullifyRe = re.compile(r" *nullify *\(", re.IGNORECASE)
    nullify2Re = re.compile(
        r"(?P<nullif> *nullify *\()(?P<vars>[^()!&]+)\)", re.IGNORECASE)

    for i in range(len(strings) - 1, -1, -1):
        line = strings[i]
        comments = []
        if nullifyRe.match(line) and findWord(var, line) != -1:
            core = ""
            comments = []
            for l in line.splitlines():
                pos = l.find("&")
                pos2 = l.find("!")
                if pos == -1:
                    if pos2 == -1:
                        core += l
                    else:
                        core += l[:pos2]
                        comments.append(l[pos2:] + "\n")
                else:
                    core += l[:pos]
                    if pos2 != -1:
                        comments.append(l[pos2:] + "\n")
            m = nullify2Re.match(core)
            if not m:
                raise SyntaxError("could not match nullify to " + repr(core) +
                                  "in" + repr(line))
            allVars = []
            vars = m.group("vars")
            v = list(map(str.strip, vars.split(",")))
            removedNow = 0
            for j in range(len(v) - 1, -1, -1):
                if findWord(var, v[j].lower()) != -1:
                    del v[j]
                    removedNow = 1
            if removedNow:
                if len(v) == 0:
                    if not comments:
                        del strings[i]
                    else:
                        strings[i] = "".join(comments)
                else:
                    for j in range(len(v) - 1):
                        v[j] += ", "
                    v[-1] += ")"
                    newS = StringIO()
                    v.insert(0, m.group("nullif"))
                    writeInCols(v, len(v[0]) -
                                len(v[0].lstrip()) + 5, 77, 0, newS)
                    newS.write("\n")
                    if comments:
                        for c in comments:
                            newS.write(c)
                    strings[i] = newS.getvalue()
                removed += 1
    return removed


def parseUse(inFile):
    """Parses the use statements in inFile
    The parsing stops at the first non use statement.
    Returns something like:
    ([{'module':'module1','only':['el1','el2=>el3']},...],
     '! comment1\\n!comment2...\\n',
     'last line (the line that stopped the parsing)')
    """
    lineNr = 0
    preComments = []
    modules = []
    origLines = []
    commonUses = ""
    stream = InputStream(inFile)
    while 1:
        (jline, comment_list, lines) = stream.nextFortranLine()
        comments = '\n'.join(_ for _ in comment_list if _)
        lineNr = lineNr + len(lines)
        if not lines:
            break
        origLines.append("".join(lines))
        # parse use
        m = USE_PARSE_RE.match(jline)
        if m:
            useAtt = {'module': m.group('module'), 'comments': [], 'intrinsic': False}

            if m.group('intrinsic'):
                useAtt['intrinsic'] = True
            if m.group('only'):
                useAtt['only'] = list(map(str.strip,
                                      str.split(m.group('imports'), ',')))
            else:
                useAtt['renames'] = list(map(str.strip,
                                         str.split(m.group('imports'), ',')))
                if useAtt['renames'] == [""]:
                    del useAtt['renames']
            if comments:
                useAtt['comments'].append(comments)
            # add use to modules
            modules.append(useAtt)
        elif jline and not jline.isspace():
            break
        else:
            if comments and commonUsesRe.match(comments):
                commonUses += "".join(lines)
            elif len(modules) == 0:
                preComments.append(("".join(lines)))
            elif comments:
                modules[-1]['comments'].append(comments)

    return {'modules': modules, 'preComments': preComments, 'commonUses': commonUses,
            'postLine': "".join(lines), 'origLines': origLines[:-1]}


def normalizeModules(modules):
    """Sorts the modules and their export and removes duplicates.
    renames aren't sorted correctly"""
    # orders modules
    modules.sort(key=lambda x: x['module'])
    for i in range(len(modules) - 1, 0, -1):
        if modules[i]['module'].lower() == modules[i - 1]['module'].lower():
            if not ('only' in modules[i - 1].keys() and
                    'only' in modules[i].keys()):
                raise SyntaxError('rejoining of module ' +
                                  str(modules[i]['module']) +
                                  ' failed as at least one of the use is not a use ...,only:')
            modules[i - 1]['only'].extend(modules[i]['only'])
            del modules[i]
    # orders imports
    for m in modules:
        if 'only' in m.keys():
            m['only'].sort()
            for i in range(len(m['only']) - 1, 0, -1):
                if m['only'][i - 1].lower() == m['only'][i].lower():
                    del m['only'][i]


def writeUses(modules, outFile):
    """Writes the use declaration using a long or short form depending on how
    many only statements there are"""
    maxoffset = 0
    lint = 0
    for m in modules:
        maxoffset = max(maxoffset,len(m['module']))
        if m['intrinsic']:
            lint = 14


    for m in modules:
        if 'only' in m.keys() and len(m['only']) > 11:
            writeUseShort(m, outFile,maxoffset+11+lint)
        else:
            writeUseLong(m, outFile,maxoffset+11+lint)


def writeUseLong(m, outFile,maxoffset):
    """Writes a use declaration in a nicer, but longer way"""
    lint = 0
    if 'only' in m.keys():
        if m['intrinsic']:
            outFile.write(INDENT_SIZE * ' ' + "Use, Intrinsic :: " + m['module'] + "," +
                      str.rjust('Only: ', maxoffset - 18 - len(m['module'])))
        else:
            outFile.write(INDENT_SIZE * ' ' + "Use " + m['module'] + "," +
                      str.rjust('Only: ', maxoffset - 4 - len(m['module'])))
        if m['only']:
            outFile.write(m['only'][0])
        for i in range(1, len(m['only'])):
            outFile.write(",&\n" + str.ljust("", maxoffset + lint + 1 +
                                                INDENT_SIZE) + m['only'][i])
    else:
        if m['intrinsic']:
            outFile.write(INDENT_SIZE * ' ' + "Use, Intrinsic :: " + m['module'])
        else:
            outFile.write(INDENT_SIZE * ' ' + "Use " + m['module'])
        if 'renames' in m.keys() and m['renames']:
            outFile.write("," + str.ljust("", maxoffset - 4) +
                          m['renames'][0])
            for i in range(1, len(m['renames'])):
                outFile.write(",&\n" + str.ljust("", maxoffset + 1 +
                                                    INDENT_SIZE) + m['renames'][i])
    if m['comments']:
        outFile.write("\n")
        outFile.write('\n'.join(m['comments']))
    outFile.write("\n")


def writeUseShort(m, file,maxoffset):
    """Writes a use declaration in a compact way"""
    uLine = []
    if 'only' in m.keys():
        if m['intrinsic']:
            file.write(INDENT_SIZE * ' ' + "Use, Intrinsic :: " + m['module'] + "," +
                   str.rjust('Only: &\n', maxoffset - 2 - len(m['module'])))
        else:
            file.write(INDENT_SIZE * ' ' + "Use " + m['module'] + "," +
                   str.rjust('Only: &\n', maxoffset - 2 - len(m['module'])))
        for k in m['only'][:-1]:
            uLine.append(k + ", ")
        uLine.append(m['only'][-1])
        uLine[0] = " " * (maxoffset + 1 + INDENT_SIZE) + uLine[0]
    elif 'renames' in m.keys() and m['renames']:
        uLine.append(INDENT_SIZE * ' ' + "Use " + m['module'] + ", ")
        for k in m['renames'][:-1]:
            uLine.append(k + ", ")
        uLine.append(m['renames'][-1])
    else:
        if m['intrinsic']:
            uLine.append(INDENT_SIZE * ' ' + "Use, Intrinsic :: " + m['module'])
        else:
            uLine.append(INDENT_SIZE * ' ' + "Use " + m['module'])
    writeInCols(uLine, maxoffset + 1 + INDENT_SIZE, DECL_LINELENGTH, 0, file)
    if m['comments']:
        file.write("\n")
        file.write('\n'.join(m['comments']))
    file.write("\n")


def prepareImplicitUses(modules):
    """Transforms a modulesDict into an implictUses (dictionary of module names
    each containing a dictionary with the only, and the special key '_WHOLE_'
    wich is true if the whole mosule is implicitly present"""
    mods = {}
    for m in modules:
        m_name = m['module'].lower()
        if m_name not in mods.keys():
            mods[m['module']] = {'_WHOLE_': 0}
        m_att = mods[m_name]
        if 'only' in m.keys():
            for k in m['only']:
                m = localNameRe.match(k)
                if not m:
                    raise SyntaxError('could not parse use only:' + repr(k))
                impAtt = m.group('localName').lower()
                m_att[impAtt] = 1
        else:
            m_att['_WHOLE_'] = 1
    return mods


def cleanUse(modulesDict, rest, implicitUses=None):
    """Removes the unneded modules (the ones that are not used in rest)"""
    logger = logging.getLogger('prettify-logger')

    global R_USE
    exceptions = set()
    modules = modulesDict['modules']
    rest = rest.lower()
    for i in range(len(modules) - 1, -1, -1):
        m_att = {}
        m_name = modules[i]['module'].lower()
        if implicitUses and m_name in implicitUses.keys():
            m_att = implicitUses[m_name]
        if '_WHOLE_' in m_att.keys() and m_att['_WHOLE_']:
            R_USE += 1
            logger.info("removed USE of module " + m_name + "\n")
            del modules[i]
        elif 'only' in modules[i].keys():
            els = modules[i]['only']
            for j in range(len(els) - 1, -1, -1):
                m = localNameRe.match(els[j])
                n = operatormRe.match(els[j])
                if not m and not n:
                    raise SyntaxError(
                        'could not parse use only:' + repr(els[j]))
                impAtt = ''
                if m:
                   impAtt = m.group('localName').lower()
                if n:
                   impAtt = n.group('localName').capitalize()
                   exceptions.add(impAtt)

                if impAtt in m_att.keys():
                    R_USE += 1
                    logger.info("removed USE " + m_name +
                                ", only: " + repr(els[j]) + "\n")
                    del els[j]
                elif impAtt not in exceptions:
                    if findWord(impAtt, rest) == -1:
                        R_USE += 1
                        logger.info("removed USE " + m_name +
                                    ", only: " + repr(els[j]) + "\n")
                        del els[j]
            if len(modules[i]['only']) == 0:
                if modules[i]['comments']:
                    modulesDict['preComments'].extend(
                        list(map(lambda x: x + "\n", modules[i]['comments'])))
                del modules[i]


def resetModuleN(moduleName, lines):
    "resets the moduleN variable to the module name in the lines lines"
    moduleNRe = re.compile(r".*:: *moduleN *= *(['\"])[a-zA-Z_0-9]+\1",
                           flags=re.IGNORECASE)
    for i in range(len(lines)):
        lines[i] = moduleNRe.sub(
            " " * INDENT_SIZE + "CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = '" +
            moduleName + "'",
            lines[i])


def rewriteFortranFile(inFile, outFile, indent, decl_linelength, orig_filename=None):
    """rewrites the use statements and declarations of inFile to outFile.
    It sorts them and removes the repetitions."""
    import os.path
    logger = logging.getLogger('prettify-logger')

    global INDENT_SIZE
    global DECL_LINELENGTH
    INDENT_SIZE = indent

    DECL_LINELENGTH = decl_linelength

    moduleRe = re.compile(r" *(?:module|program) +(?P<moduleName>[a-zA-Z_][a-zA-Z_0-9]*) *(?:!.*)?$",
                          flags=re.IGNORECASE)
    coreLines = []
    while 1:
        line = inFile.readline()
        if not line:
            break
        if line[0] == '#':
            coreLines.append(line)
        outFile.write(line)
        m = moduleRe.match(line)
        if m:
            if not orig_filename:
                orig_filename = inFile.name
            fn = os.path.basename(orig_filename).rsplit(".", 1)[0]
            break

    try:
        modulesDict = parseUse(inFile)
        routines = []
        coreLines.append(modulesDict['postLine'])
        routine = parseRoutine(inFile)
        coreLines.extend(routine['preRoutine'])
        if m:
            resetModuleN(m.group('moduleName'), routine['preRoutine'])
        routines.append(routine)
        while routine['kind']:
            routine = parseRoutine(inFile)
            routines.append(routine)
        for routine in routines:
            cleanDeclarations(routine)  # in-place modification of 'routine'
            coreLines.extend(routine['declarations'])
            coreLines.extend(routine['strippedCore'])
        rest = "".join(coreLines)
        nonStPrep = 0
        for line in modulesDict['origLines']:
            if (re.search('^#', line) and not commonUsesRe.match(line)):
                logger.debug('noMatch ' + repr(line) +
                             '\n')  # what does it mean?
                nonStPrep = 1
        if nonStPrep:
            logger.debug(
                "use statements contains preprocessor directives, not cleaning\n")
            outFile.writelines(modulesDict['origLines'])
        else:
            implicitUses = None
            if modulesDict['commonUses']:
                try:
                    inc_fn = commonUsesRe.match(
                        modulesDict['commonUses']).group(1)
                    inc_absfn = os.path.join(
                        os.path.dirname(orig_filename), inc_fn)
                    f = open(inc_absfn, 'r')
                    implicitUsesRaw = parseUse(f)
                    f.close()
                    implicitUses = prepareImplicitUses(
                        implicitUsesRaw['modules'])
                except:
                    logger.critical(
                        "ERROR trying to parse use statements contained in common uses precompiler file " + inc_absfn + '\n')
                    raise
            cleanUse(modulesDict, rest,
                     implicitUses=implicitUses)
            normalizeModules(modulesDict['modules'])
            outFile.writelines(modulesDict['preComments'])
            writeUses(modulesDict['modules'], outFile)
            outFile.write(modulesDict['commonUses'])
            if modulesDict['modules']:
                outFile.write('\n')
        outFile.write(modulesDict['postLine'])
        for routine in routines:
            writeRoutine(routine, outFile)
    except:
        import traceback
        logger.critical('-' * 60 + "\n")
        traceback.print_exc(file=sys.stderr)
        logger.critical('-' * 60 + "\n")
        logger.critical("Processing file '" + orig_filename + "'\n")
        raise

# EOF
