import sys
# from ase import Atoms

class Log:
    """logファイルを与えた場合それをLogオブジェクトとして格納する.

    Paramaters:

    eqfile: str or Path
        | XXX_EQ_list.logファイルのパス.
    tsfile: str or Path
        | XXX_EQ_list.logファイルのパス.
    ptfile: str or Path
        | XXX_EQ_list.logファイルのパス.
        
    Examples:
    
        >>> log1 = log("XXX_EQ_list.log",ts="XXX_TS_list.log", pt="XXX_PT_list.log")
        >>> log2 = log("YYY_EQ_list.log",ts="YYY_TS_list.log", pt="YYY_PT_list.log")
        >>> marge_log = log1 + log2 # マージ
        >>> marge_log.write_log("margedEQ.log", "margedTS.log", "margedPT.log")
    """
    def __init__(self, eqfile=None, tsfile=None, ptfile=None):
        """logファイルを与えた場合それをLogオブジェクトとして格納する.

        Paramaters:

        eqfile: str or Path
            XXX_EQ_list.logファイルのパス.
        tsfile: str or Path
            XXX_EQ_list.logファイルのパス.
        ptfile: str or Path
            XXX_EQ_list.logファイルのパス.
        """
        self.eqlog = []
        self.tslog = []
        self.ptlog = []

        if eqfile:
            self.eqlog = read_log(eqfile)
        if tsfile:
            self.tslog = read_log(tsfile)
        if ptfile:
            self.ptlog = read_log(ptfile)

    def __add__(self, other):
        """Logオブジェクト同士の和を以下のように定義.
        
        | EQ,TS,PTそれぞれのリストを結合し, ナンバーを連番にする.
        | それに対応したTS,PTのコネクションのナンバーを再編成する.
        
        Returns:
            Log: 単一のLogオブジェクト

        """
        if self.eqlog == [] or other.eqlog == []:
            print("ERROR: EQlist(s) not found!")
            sys.exit()
        l = len(self.eqlog) - 1
        o = len(other.eqlog)
        newlist_eq = []
        newlist_ts = []
        newlist_pt = []

        for i in range(1,o):
            n = int(other.eqlog[i].split()[4][:-1])
            newlist_eq.append(other.eqlog[i].replace(str(n), str(n + l), 1))
        self.eqlog[-1] = self.eqlog[-1].strip()
        neweq = self.eqlog + newlist_eq

        if not self.tslog == [] and not other.tslog == []:
            tsl = len(self.tslog) - 1
            tso = len(other.tslog)
            if tsl > 2 and tso > 2:
                for i in range(1,tso):
                    n = int(other.tslog[i].split()[4][:-1])
                    other.tslog[i] = other.tslog[i].replace(str(n), str(n + tsl), 1)
                    n_ini = other.tslog[i].split()[-3]
                    n_fin = other.tslog[i].split()[-1]
                    ini_fin = f"{n_ini} - {n_fin}"
                    i_ini = int(n_ini) + l
                    i_fin = int(n_fin) + l
                    new_ini_fin = f"{str(i_ini)} - {str(i_fin)}"
                    newlist_ts.append(other.tslog[i].replace(ini_fin, new_ini_fin, 1))

            self.tslog[-1] = self.tslog[-1].strip()
            newts = self.tslog + newlist_ts
        elif self.tslog == []:
            tso = len(other.tslog)
            if tso > 2:
                for i in range(1,tso):
                    n = int(other.tslog[i].split()[4][:-1])
                    n_ini = other.tslog[i].split()[-3]
                    n_fin = other.tslog[i].split()[-1]
                    ini_fin = f"{n_ini} - {n_fin}"
                    i_ini = int(n_ini) + l
                    i_fin = int(n_fin) + l
                    new_ini_fin = f"{str(i_ini)} - {str(i_fin)}"
                    newlist_ts.append(other.tslog[i].replace(ini_fin, new_ini_fin, 1))
            newts = newlist_ts
        elif other.tslog == []:
            newts = self.tslog
        else:
            newts = []

        if not self.ptlog == [] and not other.ptlog == []:
            ptl = len(self.ptlog) - 1
            pto = len(other.ptlog)
            if ptl > 2 and pto > 2:
                for i in range(1,pto):
                    n = int(other.ptlog[i].split()[4][:-1])
                    other.ptlog[i] = other.ptlog[i].replace(str(n), str(n + ptl), 1)
                    n_ini = other.ptlog[i].split()[-3]
                    n_fin = other.ptlog[i].split()[-1]
                    ini_fin = f"{n_ini} - {n_fin}"
                    i_ini = int(n_ini) + l
                    i_fin = int(n_fin) + l
                    new_ini_fin = f"{str(i_ini)} - {str(i_fin)}"
                    newlist_pt.append(other.ptlog[i].replace(ini_fin, new_ini_fin, 1))

            self.ptlog[-1] = self.ptlog[-1].strip()
            newpt = self.ptlog + newlist_pt
        elif self.ptlog == []:
            pto = len(other.ptlog)
            if pto > 2:
                for i in range(1,pto):
                    n = int(other.ptlog[i].split()[4][:-1])
                    n_ini = other.ptlog[i].split()[-3]
                    n_fin = other.ptlog[i].split()[-1]
                    ini_fin = f"{n_ini} - {n_fin}"
                    i_ini = int(n_ini) + l
                    i_fin = int(n_fin) + l
                    new_ini_fin = f"{str(i_ini)} - {str(i_fin)}"
                    newlist_pt.append(other.ptlog[i].replace(ini_fin, new_ini_fin, 1))
            newpt = newlist_pt
        elif other.ptlog == []:
            newpt = self.ptlog
        else:
            newpt = []

        newobj = Log()
        newobj.eqlog = neweq
        newobj.tslog = newts
        newobj.ptlog = newpt
        return newobj
    
    def __iadd__(self,other):
        return self + other

    def write_log(self, eqfile="margedEQ.log", tsfile="margedTS.log", ptfile="margedPT.log"):
        """Logオブジェクトからlogファイルを書き出す.

        Paramaters:

        eqfile: str or Path
            | 保存するXXX_EQ_list.logファイルのパス.
        tsfile: str or Path
            | 保存するXXX_EQ_list.logファイルのパス.
        ptfile: str or Path
            | 保存するXXX_EQ_list.logファイルのパス.
        """
        if not eqfile == None and not self.eqlog[-2] == "List of Equilibrium Structures":
            with open(eqfile, mode='w') as f:
                for eq in self.eqlog:
                    if eq == self.eqlog[-1]:
                        f.writelines(eq)
                    else:
                        f.writelines(eq +  "\n\n")
                f.close()
        if not tsfile == None and not self.tslog == None:
            if not self.tslog[-2] == "List of Transition Structures":
                with open(tsfile, mode='w') as f:
                    for ts in self.tslog:
                        if ts == self.tslog[-1]:
                            f.writelines(ts)
                        else:
                            f.writelines(ts +  "\n\n")
                    f.close()
        if not ptfile == None and not self.ptlog == None:
            if not self.ptlog[-2] == "List of Path Top (Approximate TS) Structures":
                with open(ptfile, mode='w') as f:
                    for pt in self.ptlog:
                        if pt == self.ptlog[-1]:
                            f.writelines(pt)
                        else:
                            f.writelines(pt +  "\n\n")
                    f.close()
        
    # def eqatoms(self):
    #     l = len(self.eqlog)
    #     a_lists = []
    #     for i in range(1, l):
    #         a_list = []
    #         s = self.eqlog[i].split('\n')
    #         for j in range(1, 10000):
    #             if s[j][:6] == "Energy":
    #                 break
    #             else:
    #                 a_list.append(s[j])
            
    #         a_lists.append(a_list)
        
    #     return a_lists

def read_log(logfile):
    """logファイルを各EQ(or TS or PT)ごとに切り分けてリストで返す.
    
    Paramater:
    
    logfile: str or Path
        logファイルのパス
        
    Return:
        list of str: 各EQ(or TS or PT)のリスト"""
    log_list = open(logfile)
    r = log_list.read()
    r_list = r.split("\n\n")
    log_list.close()
    return r_list