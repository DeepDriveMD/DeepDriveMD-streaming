import os, sys, h5py
import adios2
import numpy as np
import mytimer
from aggregator_reader import *

# from keras.optimizers import RMSprop

from cvae.vae_conv import conv_variational_autoencoder
# sys.path.append('/home/hm0/Research/molecules/molecules_git/build/lib')
# from molecules.ml.unsupervised import VAE
# from molecules.ml.unsupervised import EncoderConvolution2D
# from molecules.ml.unsupervised import DecoderConvolution2D
# from molecules.ml.unsupervised.callbacks import EmbeddingCallback 

# def CVAE(input_shape, hyper_dim=3): 
#     optimizer = RMSprop(lr=0.001, rho=0.9, epsilon=1e-08, decay=0.0)

#     encoder = EncoderConvolution2D(input_shape=input_shape)

#     encoder._get_final_conv_params()
#     num_conv_params = encoder.total_conv_params
#     encode_conv_shape = encoder.final_conv_shape

#     decoder = DecoderConvolution2D(output_shape=input_shape,
#                                    enc_conv_params=num_conv_params,
#                                    enc_conv_shape=encode_conv_shape)

#     cvae = VAE(input_shape=input_shape,
#                latent_dim=hyper_dim,
#                encoder=encoder,
#                decoder=decoder,
#                optimizer=optimizer) 
#     return cvae 

def CVAE(input_shape, latent_dim=3): 
    image_size = input_shape[:-1]
    channels = input_shape[-1]
    conv_layers = 4
    feature_maps = [64,64,64,64]
    filter_shapes = [(3,3),(3,3),(3,3),(3,3)]
    strides = [(1,1),(2,2),(1,1),(1,1)]
    dense_layers = 1
    dense_neurons = [128]
    dense_dropouts = [0.25]

    feature_maps = feature_maps[0:conv_layers];
    filter_shapes = filter_shapes[0:conv_layers];
    strides = strides[0:conv_layers];
    autoencoder = conv_variational_autoencoder(image_size,channels,conv_layers,feature_maps,
               filter_shapes,strides,dense_layers,dense_neurons,dense_dropouts,latent_dim); 
    autoencoder.model.summary()
    return autoencoder

def run_cvae(gpu_id, cm_files, hyper_dim=3, epochs=20, nsamples=100, init_cvae=None): 
    # read contact map from h5 file 

    inputs = []

    mytimer.mytime_label("cvae_read", 1)

    ADIOS_XML=os.path.dirname(cm_files[0]) + "/adios.xml"

    streams = STREAMS(cm_files, lastN=nsamples, config=ADIOS_XML)
    cm_data_input = streams.next_cm()

    '''
    for cm_file in cm_files:
        print(f"cm_file={cm_file}")
        sys.stdout.flush()
        with adios2.open(cm_file, "r") as fr:
            n = fr.steps()
            print(f"n={n}"); sys.stdout.flush()
            vars = fr.available_variables()
            name = 'contact_map'
            shape = list(map(int, vars[name]['Shape'].split(",")))
            zs = list(np.zeros(len(shape), dtype='int'))
            cm_data_input = fr.read(name, zs, shape, 0, n)
    
        print(f"bp.shape = {cm_data_input.shape}")

        cm_data_input = cm_data_input[-nsamples:]
        np.random.shuffle(cm_data_input)

        inputs.append(cm_data_input)
    cm_data_input = np.concatenate(inputs)
    '''
    np.random.shuffle(cm_data_input)

    mytimer.mytime_label("cvae_read", -1)

    # splitting data into train and validation
    train_val_split = int(0.8 * len(cm_data_input))
    cm_data_train, cm_data_val = cm_data_input[:train_val_split], cm_data_input[train_val_split:] 

    input_shape = cm_data_train.shape

    # os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"
    # os.environ["CUDA_VISIBLE_DEVICES"]=str(gpu_id) 
    

    print(f"CUDA_DEVICE_ORDER = ", os.getenv("CUDA_DEVICE_ORDER"))
    print(f"CUDA_VISIBLE_DEVICES = ", os.getenv("CUDA_VISIBLE_DEVICES"))
    print(f"gpu_id = {gpu_id}")


    if(init_cvae==None):
        cvae = CVAE(input_shape[1:], hyper_dim) 
    else:
        cvae = init_cvae
    
    #     callback = EmbeddingCallback(cm_data_train, cvae)
    cvae.train(cm_data_train, validation_data=cm_data_val, batch_size = min(1024, int(input_shape[0]/100)), epochs=epochs) 
    #cvae.train(cm_data_train, validation_data=cm_data_val, batch_size = 128, epochs=epochs) 
    
    return cvae 


