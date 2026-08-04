import re

input_conf = "24x-heatmap-circos-sumai3.conf"
output_conf = "60x-heatmap-circos-sumai3.conf"

scale = 3.0   # 2000p -> 6000p


with open(input_conf, "r") as f:
    text = f.read()


# -----------------------------
# Scale all pixel values (p)
# -----------------------------
def scale_p(match):
    value = float(match.group(1))
    new_value = value * scale

    if new_value.is_integer():
        return f"{int(new_value)}p"
    else:
        return f"{new_value:.2f}p"


text = re.sub(
    r'(\d+(?:\.\d+)?)p',
    scale_p,
    text
)


# -----------------------------
# Scale label_size without p
# -----------------------------
def scale_label(match):
    value = float(match.group(1))
    new_value = value * scale

    if new_value.is_integer():
        return f"label_size       = {int(new_value)}"
    else:
        return f"label_size       = {new_value:.2f}"


text = re.sub(
    r'label_size\s+=\s+(\d+(?:\.\d+)?)\s*$',
    scale_label,
    text,
    flags=re.MULTILINE
)


# -----------------------------
# Rename output
# -----------------------------
text = text.replace(
    "file  = 24circos-sumai3.png",
    "file  = 60circos-sumai3.png"
)


# -----------------------------
# Save config
# -----------------------------
with open(output_conf, "w") as f:
    f.write(text)


print("Saved:", output_conf)
print("Scaled 2000p -> 6000p")
