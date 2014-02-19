
# Tomorrow palette
# https://github.com/chriskempson/tomorrow-theme
systemic.theme.tomorrow <- c(
    'black',
    '#4271ae',
    '#c82829',
    '#718c00',
    '#8959a8',
    '#3e999f',
    '#f5871f',
    '#eab700'
)

systemic.theme.tomorrow.face <- adjustcolor(systemic.theme.tomorrow, offset=c(0.2, 0.2, 0.2, 0.2))

# Default R palette
systemic.theme.r <- palette("default")
systemic.theme.r.face <- palette("default")

systemic.palette <- systemic.theme.tomorrow
systemic.palette.face <- systemic.theme.tomorrow.face
