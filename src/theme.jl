black = :gray10
green = :seagreen1
green_dimmed = :seagreen
green_lessdimmed = :mediumseagreen
colormap = :hawaii
font = "FreeMono Bold"

kstheme = Theme(;
    backgroundcolor=black,
    textcolor=green,
    linecolor=:white,
    fonts=(; regular=font, bold=font),
    fontsize=14,
    Axis=(titlesize=16, backgroundcolor=:transparent,
          bottomspinecolor=green, topspinecolor=green, leftspinecolor=green,
          rightspinecolor=green, xgridcolor=RGBAf(1, 1, 1, 0.16),
          ygridcolor=RGBAf(1, 1, 1, 0.16), xtickcolor=green, ytickcolor=green),
    Legend=(framecolor=green, bgcolor=black),
    Slider=(color_active=green, color_active_dimmed=green_dimmed,
            color_inactive=green_dimmed),
    Colorbar=(tickcolor=green, spinecolor=green, topspinecolor=green,
              bottomspinecolor=green, leftspinecolor=green, rightspinecolor=green,
              colormap=colormap),
    Toggle=(buttoncolor=green, toggleduration=0.1,),
    # TODO: breaks the toggle? framecolor_active=:black, framecolor_inactive=:white), 
    Textbox=(bordercolor=green, bordercolor_focused=green,
             bordercolor_hover=green_lessdimmed),
    Button=(buttoncolor=green_dimmed, buttoncolor_active=green_lessdimmed),
)
