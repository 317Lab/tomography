#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#
# SPDX-License-Identifier: GPL-3.0
#
# GNU Radio Python Flow Graph
# Title: 317 IQ Pre-processing
# GNU Radio version: 3.10.9.2

from gnuradio import blocks
import pmt
from gnuradio import filter
from gnuradio.filter import firdes
from gnuradio import gr
from gnuradio.fft import window
import sys
import signal
from argparse import ArgumentParser
from gnuradio.eng_arg import eng_float, intx
from gnuradio import eng_notation




class decsink(gr.top_block):

    def __init__(self, sink_name='0', source_name='0', startind=0):
        gr.top_block.__init__(self, "317 IQ Pre-processing", catch_exceptions=True)

        ##################################################
        # Parameters
        ##################################################
        self.sink_name = sink_name
        self.source_name = source_name
        self.startind = startind

        ##################################################
        # Variables
        ##################################################
        self.samp_rate = samp_rate = 50000
        self.dec = dec = 16

        ##################################################
        # Blocks
        ##################################################

        self.filter_fft_low_pass_filter_0 = filter.fft_filter_ccc(dec, firdes.low_pass(1, samp_rate, (samp_rate/(dec*2)), (samp_rate/(dec*2*16)), window.WIN_HANN, 6.76), 1)
        self.blocks_file_source_0_0 = blocks.file_source(gr.sizeof_gr_complex*1, source_name, False, startind, 0)
        self.blocks_file_source_0_0.set_begin_tag(pmt.PMT_NIL)
        self.blocks_file_sink_0 = blocks.file_sink(gr.sizeof_gr_complex*1, sink_name, False)
        self.blocks_file_sink_0.set_unbuffered(False)


        ##################################################
        # Connections
        ##################################################
        self.connect((self.blocks_file_source_0_0, 0), (self.filter_fft_low_pass_filter_0, 0))
        self.connect((self.filter_fft_low_pass_filter_0, 0), (self.blocks_file_sink_0, 0))


    def get_sink_name(self):
        return self.sink_name

    def set_sink_name(self, sink_name):
        self.sink_name = sink_name
        self.blocks_file_sink_0.open(self.sink_name)

    def get_source_name(self):
        return self.source_name

    def set_source_name(self, source_name):
        self.source_name = source_name
        self.blocks_file_source_0_0.open(self.source_name, False)

    def get_startind(self):
        return self.startind

    def set_startind(self, startind):
        self.startind = startind

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.filter_fft_low_pass_filter_0.set_taps(firdes.low_pass(1, self.samp_rate, (self.samp_rate/(self.dec*2)), (self.samp_rate/(self.dec*2*16)), window.WIN_HANN, 6.76))

    def get_dec(self):
        return self.dec

    def set_dec(self, dec):
        self.dec = dec
        self.filter_fft_low_pass_filter_0.set_taps(firdes.low_pass(1, self.samp_rate, (self.samp_rate/(self.dec*2)), (self.samp_rate/(self.dec*2*16)), window.WIN_HANN, 6.76))



def argument_parser():
    parser = ArgumentParser()
    parser.add_argument(
        "--sink-name", dest="sink_name", type=str, default='0',
        help="Set default_processed.bin [default=%(default)r]")
    parser.add_argument(
        "--source-name", dest="source_name", type=str, default='0',
        help="Set default.bin [default=%(default)r]")
    parser.add_argument(
        "--startind", dest="startind", type=intx, default=0,
        help="Set startind [default=%(default)r]")
    return parser


def main(top_block_cls=decsink, options=None):
    if options is None:
        options = argument_parser().parse_args()
    tb = top_block_cls(sink_name=options.sink_name, source_name=options.source_name, startind=options.startind)

    def sig_handler(sig=None, frame=None):
        tb.stop()
        tb.wait()

        sys.exit(0)

    signal.signal(signal.SIGINT, sig_handler)
    signal.signal(signal.SIGTERM, sig_handler)

    tb.start()

    tb.wait()


if __name__ == '__main__':
    main()
