#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
out_plot_prefix <- args[2]
out_dir <- args[3]
format <- args[4]  # "png" ou "pdf"
mode <- ifelse(length(args) >= 5, args[5], "compact")

# Créer le dossier de sortie si besoin
if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
}

# Packages
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(readxl))
suppressMessages(library(RColorBrewer))
suppressMessages(library(stringr))


# Fonction pour enregistrer un plot
save_plot <- function(plot, filename, format, width, height) {
    filepath <- file.path(out_dir, paste0(filename, ".", format))
    
    # Utiliser dpi uniquement pour PNG
    if (format == "png") {
        ggsave(
            filepath,
            plot = plot,
            width = width,
            height = height,
            dpi = 300,  # résolution de 300 pour PNG
            bg = "white"  # Fond blanc garanti
        )
    } else {
        ggsave(
            filepath,
            plot = plot,
            width = width,
            height = height,
            bg = "white"  # Fond blanc garanti, pas de dpi pour PDF
        )
    }
}



# Lire le fichier
data <- read_excel(file)

# Créer colonnes Type / Lineage / Sub-lineage selon le mode
if (mode == "compact") {
    data$`PAVE match` <- trimws(data$`PAVE match`)
    data <- data[data$`PAVE match` != "ND" & data$`PAVE match` != "", ]

    splitted <- strsplit(as.character(data$`PAVE match`), "-")
    data_clean <- data.frame(
        Type = sapply(splitted, `[`, 1),
        Lineage = sapply(splitted, `[`, 2),
        Sub_lineage = sapply(splitted, `[`, 3),
        stringsAsFactors = FALSE
    )

} else if (mode == "detailed") {
    Type_col <- 10  # J
    Lineage_col <- 11  # K
    Sublineage_col <- 12  # L

    data_clean <- data.frame(
        Type = trimws(data[[Type_col]]),
        Lineage = trimws(data[[Lineage_col]]),
        Sub_lineage = trimws(data[[Sublineage_col]]),
        stringsAsFactors = FALSE
    )

    data_clean <- data_clean %>%
        filter(Type != "", !is.na(Type))  # garder seulement ceux avec un Type

} else {
    stop("Mode non reconnu : utilisez 'compact' ou 'detailed'.")
}

# Nettoyage supplémentaire
data_clean <- data_clean %>%
    mutate(
        Type = gsub("[^A-Za-z0-9_-]", "", Type),
        Type = replace(Type, Type == "" | Type == "ND", NA),
        Lineage = trimws(Lineage) %>% replace(. %in% c("", "ND"), NA),
        Sub_lineage = trimws(Sub_lineage) %>% replace(. %in% c("", "ND"), NA)
    ) %>%
    filter(!is.na(Type))

# ------------------------------------------
# Thème ggplot
# ------------------------------------------
custom_theme <- function(base_size = 14) {
    theme_bw(base_size = base_size) +
    theme(
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(color = "black"),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white")
    )
}

# ------------------------------------------
# Plot 1 : Distribution des Types
# ------------------------------------------
top_n <- 20
type_counts <- data_clean %>%
    count(Type, sort = TRUE) %>%
    head(top_n)

p1 <- ggplot(type_counts, aes(x = reorder(Type, -n), y = n)) +
    geom_col(fill = "steelblue", show.legend = FALSE) +
    geom_text(aes(label = n), vjust = -0.5, size = 3, color = "black") +
    labs(
        title = if (nrow(type_counts) < top_n) paste("Top", nrow(type_counts), "Types") else "Top 20 Types",
        x = "Type", 
        y = "Number of Sequences"
    ) +
    custom_theme()

save_plot(p1, paste0(out_plot_prefix, "_types"), format, 8, 5)

# ------------------------------------------
# Plot 2 : Distribution des Lineages par Type
# ------------------------------------------
# Distribution des Lineages par Type (LIMITÉ AUX TOP N TYPES)

top_types <- type_counts$Type  # on reprend les top types du premier plot

lineage_data <- data_clean %>%
    filter(Type %in% top_types, !is.na(Lineage)) %>%
    count(Type, Lineage)

p2 <- ggplot(lineage_data, aes(x = Type, y = n, fill = Lineage)) +
    geom_col(position = "dodge") +
    labs(
        title = paste0("Lineages per HPV Type (Top ", length(top_types), ")"),
        x = "Type",
        y = "Number of Sequences"
    ) +
    scale_fill_viridis_d() +
    custom_theme()

save_plot(p2, paste0(out_plot_prefix, "_lineages_per_type"), format, 10, 6)


# ------------------------------------------
# Plot 3 : Distribution des Sub-lineages par Lineage (Top 10 Types)
# ------------------------------------------
top_types_sublineage <- data_clean %>%
    filter(!is.na(Type), !is.na(Lineage), !is.na(Sub_lineage)) %>%
    count(Type, Lineage, Sub_lineage) %>%
    group_by(Type) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    top_n(10, n) %>%
    pull(Type)

sublineage_data_top <- data_clean %>%
    filter(Type %in% top_types_sublineage, !is.na(Lineage), !is.na(Sub_lineage)) %>%
    count(Type, Lineage, Sub_lineage)


# Extraire chiffre de fin (ex: A1 → 1, B2 → 2)
sublineage_data_top <- sublineage_data_top %>%
  mutate(SublineageGroup = str_extract(Sub_lineage, "\\d+"))


# Nombre unique de groupes
n_groups <- length(unique(sublineage_data_top$SublineageGroup))

# Générer palette répétée ou interpolée si besoin
if (n_groups <= 8) {
  my_colors <- brewer.pal(n_groups, "Dark2")
} else {
  my_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(n_groups)
}



p3 <- ggplot(sublineage_data_top, aes(x = Lineage, y = n, fill = SublineageGroup)) +
  geom_col(position = "dodge") +
  facet_wrap(~ Type, scales = "free_x") +
  scale_fill_manual(values = my_colors) +  # ggplot utilisera les couleurs dans l'ordre
  labs(
    title = "Top 10 Types - Sublineages per Lineage",
    x = "Lineage",
    y = "Number of Sequences"
  ) +
  custom_theme()

save_plot(p3, paste0(out_plot_prefix, "_sublineages_top10_types"), format, 14, 8)

# ------------------------------------------

