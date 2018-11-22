.libPaths('~/R_libs')

library(bookdown)
bookdown::render_book(input = 'index.Rmd', ## all files that match '.Rmd' as input
                      # output_format = 'all', ## only gitbook outbook, no pdf or ePub
                      output_format = 'bookdown::gitbook',
                      new_session = T, ## TRUE: do not merge the Rmd files into a single file
                      preview=F, ## ???
                      output_dir= '../docs', ## where to store the output vignettes
                      params =   list( ## input list of parameters
                        local.input='../data', ## input data - F for github load
                        subset.data = F, ## a reduced set of data; or F for main data
                        # subset.data = '../data/subset/sincell_reduced_80_200.RData',
                        output.data='../output', ## to store saved data
                        Rscripts = '../Rscripts', ## where to save r scripts - F for not saving
                        recalc = F ## do not recalculate the objects for second vignette
                      ),
                      clean=T
                      
)
## clean up
file.remove(list.files(pattern = '.rds'))
file.remove(list.files(pattern = '.log'))
file.remove(list.files(pattern = '.png'))
## add the necessary files to stand-alone folder
file.copy('01-Data-Integration.Rmd', to = '../vignettes-standalone/01-Data-Integration.Rmd', overwrite = T)
file.copy('02-Signature.Rmd', to = '../vignettes-standalone/02-Signature.Rmd', overwrite = T)
file.copy(list.files(pattern = '.bib'), to = '../vignettes-standalone',recursive = T, overwrite = T)
file.copy('figures', to = '../vignettes-standalone', recursive = T, overwrite = T)
