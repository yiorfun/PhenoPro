readNews <- function () 
{
    newsfile <- file.path(system.file(package = "PhenoPro"), 
        "NEWS")
    file.show(newsfile)
}
