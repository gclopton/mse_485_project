# scripts/yaml2lmpvars.py
#!/usr/bin/env python3
import sys, yaml

def flat(prefix, obj, out):
    if isinstance(obj, dict):
        for k, v in obj.items(): flat(f"{prefix}_{k}".strip("_"), v, out)
    elif isinstance(obj, (list, tuple)):
        pass
    else:
        out[prefix] = obj

def is_num(x):
    try: float(x); return True
    except: return False

if __name__ == "__main__":
    if len(sys.argv)!=3:
        sys.exit("Usage: yaml2lmpvars.py <physics.yaml> <out.varfile>")
    src, dst = sys.argv[1], sys.argv[2]
    data = yaml.safe_load(open(src))
    f = {}
    flat("", data, f)
    with open(dst, "w") as o:
        o.write(f"# generated from {src}\n")
        for k in sorted(f.keys()):
            if not k: continue
            name = "p_"+k.replace("__","_")
            v = f[k]
            if isinstance(v, bool):
                o.write(f"variable {name} equal {1 if v else 0}\n")
            elif is_num(v):
                o.write(f"variable {name} equal {v}\n")
            else:
                s = str(v).replace('"','\\"')
                o.write(f'variable {name} index "{s}"\n')

