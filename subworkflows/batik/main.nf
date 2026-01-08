/*

Copyright (c) 2025 Pierre Lindenbaum

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
The MIT License (MIT)
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
include { BATIK_DOWNLOAD                          } from '../../modules/batik/download'
include { BATIK_RASTERIZE as  BATIK_RASTERIZE_PDF } from '../../modules/batik/rasterize'
include { BATIK_RASTERIZE as  BATIK_RASTERIZE_PNG } from '../../modules/batik/rasterize'
include { BATIK_RASTERIZE as  BATIK_RASTERIZE_JPG } from '../../modules/batik/rasterize'
include { parseBoolean                            } from '../../modules/utils/functions'



workflow BATIK {
take:
	metadata
	svgs
main:
	versions = Channel.empty()
    pdf = Channel.empty()
    png = Channel.empty()
    jpg = Channel.empty()

    BATIK_DOWNLOAD( svgs.collect().map{[id:"batik"]} )
    versions = versions.mix(BATIK_DOWNLOAD.out.versions)

    

    if(metadata.with_pdf==null || parseBoolean(metadata.with_pdf)) {
        BATIK_RASTERIZE_PDF(
            BATIK_DOWNLOAD.out.rasterizer,
            svgs
            )
        versions = versions.mix(BATIK_RASTERIZE_PDF.out.versions)
        pdf = pdf.mix(BATIK_RASTERIZE_PDF.out.image)
        }
    
     if(metadata.with_png==null || parseBoolean(metadata.with_png)) {
        BATIK_RASTERIZE_PNG(
            BATIK_DOWNLOAD.out.rasterizer,
            svgs
            )
        versions = versions.mix(BATIK_RASTERIZE_PNG.out.versions)
        png = pdf.mix(BATIK_RASTERIZE_PNG.out.image)
        }
     if(metadata.with_jpg==null || parseBoolean(metadata.with_jpg)) {
        BATIK_RASTERIZE_JPG(
            BATIK_DOWNLOAD.out.rasterizer,
            svgs
            )
        versions = versions.mix(BATIK_RASTERIZE_JPG.out.versions)
        jpg = jpg.mix(BATIK_RASTERIZE_JPG.out.image)
        }
emit:
    versions
    pdf
    png
    jpg
}