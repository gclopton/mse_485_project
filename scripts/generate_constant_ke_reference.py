# Create a constant ke_reference.csv for Phase 0
import pandas as pd

# Constant ke (W/m·K) and temperature grid (K)
ke_value = 1.5
temps = list(range(300, 6301, 300))

df = pd.DataFrame({"T_K": temps, "ke_W_mK": [ke_value]*len(temps)})
path = "inputs/ttm/ke_reference.csv"
df.to_csv(path, index=False)

# Compute g0 using r0 = 4.5 nm
r0_m = 4.5e-9
g0 = ke_value/(r0_m**2)

# If successfull print message to terminal
print(f"Constant ke_reference.csv created at {path}")
print(f"g0 = {g0} W/m^3/K")
