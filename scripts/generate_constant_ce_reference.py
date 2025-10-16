# Create a constant Ce_reference.csv for Phase 0
import pandas as pd

# Choose a constant Ce (J/m^3/K) to yield reasonable electron diffusivity De = ke/Ce
Ce_value = 5.0e4  # J/m^3/K
temps = list(range(300, 6301, 300))
df = pd.DataFrame({"T_K": temps, "Ce_J_m3K": [Ce_value]*len(temps)})
path = "inputs/ttm/Ce_reference.csv"
df.to_csv(path, index=False)

# If successfull print message to terminal
print(f"Constant Ce_reference.csv created at {path}")
