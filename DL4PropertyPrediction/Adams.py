"""
The MIT License (MIT)

Copyright (c) 2015 Alec Radford

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""  

import theano
import theano.tensor as T
import numpy as np
  

def Adam(params, grads, lr=0.0002, b1=0.1, b2=0.001, e=1e-8):
    updates = []
    other_params = []

    #grads = T.grad(cost, params)
    i = theano.shared(np.float32(0.).astype(theano.config.floatX))
    other_params.append(i)
    i_t = i + 1.

    fix1 = 1. - (1. - b1)**i_t
    fix2 = 1. - (1. - b2)**i_t
    lr_t = lr * (T.sqrt(fix2) / fix1)
    for p, g in zip(params, grads):
        m = theano.shared(p.get_value() * 0.)
        v = theano.shared(p.get_value() * 0.)
	other_params.append(m)
	other_params.append(v)

        m_t = (b1 * g) + ((1. - b1) * m)
        v_t = (b2 * T.sqr(g)) + ((1. - b2) * v)
        g_t = m_t / (T.sqrt(v_t) + e)
        p_t = p - (lr_t * g_t)

        updates.append((m, m_t))
        updates.append((v, v_t))
        updates.append((p, p_t))

    updates.append((i, i_t))
    return updates, other_params

def AMSGrad(params, grads, lr=0.0002, b1=0.1, b2=0.001, e=1e-8):
    updates = []
    other_params = []

    #grads = T.grad(cost, params)
    i = theano.shared(np.float32(0.).astype(theano.config.floatX))
    other_params.append(i)
    i_t = i + 1.

    fix1 = 1. - (1. - b1)**i_t
    fix2 = 1. - (1. - b2)**i_t
    lr_t = lr * (T.sqrt(fix2) / fix1)
    for p, g in zip(params, grads):
        m = theano.shared(p.get_value() * 0.)
        v = theano.shared(p.get_value() * 0.)
        v_hat = theano.shared(p.get_value() * 0.)
	other_params.append(m)
	other_params.append(v)
	other_params.append(v_hat)

        m_t = (b1 * g) + ((1. - b1) * m)
        v_t = (b2 * T.sqr(g)) + ((1. - b2) * v)
	v_hat_t = T.maximum( v_hat, v_t)

        g_t = m_t / (T.sqrt(v_hat_t) + e)
        p_t = p - (lr_t * g_t)

        updates.append((m, m_t))
        updates.append((v, v_t))
	updates.append((v_hat, v_hat_t) )
        updates.append((p, p_t))

    updates.append((i, i_t))
    return updates, other_params

## pdecay has the same length of params. Each element in pdecay is either 0 or corresponds to params
def AdamW(params, grads, pdecay=None, l2reg=0.1, lr=0.0002, b1=0.1, b2=0.001, e=1e-8):
    if pdecay is None:
	wdecay = params
    else:
	wdecay = pdecay

    updates = []
    other_params = []

    #grads = T.grad(cost, params)
    i = theano.shared(np.float32(0.).astype(theano.config.floatX))
    other_params.append(i)
    i_t = i + 1.

    fix1 = 1. - (1. - b1)**i_t
    fix2 = 1. - (1. - b2)**i_t
    lr_t = lr * (T.sqrt(fix2) / fix1)
    for p, g, d in zip(params, grads, wdecay):

	m =  theano.shared(value=np.zeros(p.shape.eval(),dtype=theano.config.floatX), borrow=True)
	v =  theano.shared(value=np.zeros(p.shape.eval(),dtype=theano.config.floatX), borrow=True)

	other_params.append(m)
	other_params.append(v)

        m_t = (b1 * g) + ((1. - b1) * m)
        v_t = (b2 * T.sqr(g)) + ((1. - b2) * v)
        g_t = m_t / (T.sqrt(v_t) + e)
	
        p_t = p - lr_t * g_t - lr * l2reg * d 

        updates.append((m, m_t))
        updates.append((v, v_t))
        updates.append((p, p_t))

    updates.append((i, i_t))
    return updates, other_params

## pdecay has the same length of params. Each element in pdecay is either 0 or correspond to params
def AdamWAMS(params, grads, pdecay=None, l2reg=0.1, lr=0.0002, b1=0.1, b2=0.001, e=1e-8):
    if pdecay is None:
	wdecay = params
    else:
	wdecay = pdecay

    updates = []
    other_params = []

    #grads = T.grad(cost, params)
    i = theano.shared(np.float32(0.).astype(theano.config.floatX))
    other_params.append(i)
    i_t = i + 1.

    fix1 = 1. - (1. - b1)**i_t
    fix2 = 1. - (1. - b2)**i_t
    lr_t = lr * (T.sqrt(fix2) / fix1)
    for p, g, d in zip(params, grads, wdecay):

        m =  theano.shared(value=np.zeros(p.shape.eval(),dtype=theano.config.floatX), borrow=True)
        v =  theano.shared(value=np.zeros(p.shape.eval(),dtype=theano.config.floatX), borrow=True)
        v_hat =  theano.shared(value=np.zeros(p.shape.eval(),dtype=theano.config.floatX), borrow=True)

        other_params.append(m)
        other_params.append(v)
        other_params.append(v_hat)

        m_t = (b1 * g) + ((1. - b1) * m)
        v_t = (b2 * T.sqr(g)) + ((1. - b2) * v)
	v_hat_t = T.maximum( v_hat, v_t)

        g_t = m_t / (T.sqrt(v_hat_t) + e)
        p_t = p - lr_t * g_t - lr * l2reg * d

        updates.append((m, m_t))
        updates.append((v, v_t))
	updates.append((v_hat, v_hat_t) )
        updates.append((p, p_t))

    updates.append((i, i_t))
    return updates, other_params

