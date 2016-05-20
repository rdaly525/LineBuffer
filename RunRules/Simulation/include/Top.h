/**\file */
#ifndef SLIC_DECLARATIONS_Top_H
#define SLIC_DECLARATIONS_Top_H
#include "MaxSLiCInterface.h"
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#define Top_nxMax (4096)
#define Top_pad (1)
#define Top_stride (1)
#define Top_k (3)
#define Top_PCIE_ALIGNMENT (16)


/*----------------------------------------------------------------------------*/
/*---------------------------- Interface readLMem ----------------------------*/
/*----------------------------------------------------------------------------*/



/**
 * \brief Basic static function for the interface 'readLMem'.
 * 
 * \param [in] param_size Interface Parameter "size".
 * \param [in] param_start Interface Parameter "start".
 * \param [out] outstream_tocpu The stream should be of size param_size bytes.
 */
void Top_readLMem(
	uint32_t param_size,
	uint32_t param_start,
	float *outstream_tocpu);

/**
 * \brief Basic static non-blocking function for the interface 'readLMem'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] param_size Interface Parameter "size".
 * \param [in] param_start Interface Parameter "start".
 * \param [out] outstream_tocpu The stream should be of size param_size bytes.
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *Top_readLMem_nonblock(
	uint32_t param_size,
	uint32_t param_start,
	float *outstream_tocpu);

/**
 * \brief Advanced static interface, structure for the engine interface 'readLMem'
 * 
 */
typedef struct { 
	uint32_t param_size; /**<  [in] Interface Parameter "size". */
	uint32_t param_start; /**<  [in] Interface Parameter "start". */
	float *outstream_tocpu; /**<  [out] The stream should be of size param_size bytes. */
} Top_readLMem_actions_t;

/**
 * \brief Advanced static function for the interface 'readLMem'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void Top_readLMem_run(
	max_engine_t *engine,
	Top_readLMem_actions_t *interface_actions);

/**
 * \brief Advanced static non-blocking function for the interface 'readLMem'.
 *
 * Schedule the actions to run on the engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in] interface_actions Actions to be executed.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Top_readLMem_run_nonblock(
	max_engine_t *engine,
	Top_readLMem_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'readLMem'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void Top_readLMem_run_group(max_group_t *group, Top_readLMem_actions_t *interface_actions);

/**
 * \brief Group run advanced static non-blocking function for the interface 'readLMem'.
 * 
 *
 * Schedule the actions to run on the first device available in the group and return immediately.
 * The status of the run must be checked with ::max_wait. 
 * Note that use of ::max_nowait is prohibited with non-blocking running on groups:
 * see the ::max_run_group_nonblock documentation for more explanation.
 *
 * \param [in] group Group to use.
 * \param [in] interface_actions Actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Top_readLMem_run_group_nonblock(max_group_t *group, Top_readLMem_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'readLMem'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void Top_readLMem_run_array(max_engarray_t *engarray, Top_readLMem_actions_t *interface_actions[]);

/**
 * \brief Array run advanced static non-blocking function for the interface 'readLMem'.
 * 
 *
 * Schedule to run the array of actions on the array of engines, and return immediately.
 * The length of interface_actions must match the size of engarray.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * \param [in] engarray The array of devices to use.
 * \param [in] interface_actions The array of actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Top_readLMem_run_array_nonblock(max_engarray_t *engarray, Top_readLMem_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* Top_readLMem_convert(max_file_t *maxfile, Top_readLMem_actions_t *interface_actions);



/*----------------------------------------------------------------------------*/
/*--------------------------- Interface writeLMem ----------------------------*/
/*----------------------------------------------------------------------------*/



/**
 * \brief Basic static function for the interface 'writeLMem'.
 * 
 * \param [in] param_size Interface Parameter "size".
 * \param [in] param_start Interface Parameter "start".
 * \param [in] instream_fromcpu The stream should be of size param_size bytes.
 */
void Top_writeLMem(
	uint32_t param_size,
	uint32_t param_start,
	const float *instream_fromcpu);

/**
 * \brief Basic static non-blocking function for the interface 'writeLMem'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] param_size Interface Parameter "size".
 * \param [in] param_start Interface Parameter "start".
 * \param [in] instream_fromcpu The stream should be of size param_size bytes.
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *Top_writeLMem_nonblock(
	uint32_t param_size,
	uint32_t param_start,
	const float *instream_fromcpu);

/**
 * \brief Advanced static interface, structure for the engine interface 'writeLMem'
 * 
 */
typedef struct { 
	uint32_t param_size; /**<  [in] Interface Parameter "size". */
	uint32_t param_start; /**<  [in] Interface Parameter "start". */
	const float *instream_fromcpu; /**<  [in] The stream should be of size param_size bytes. */
} Top_writeLMem_actions_t;

/**
 * \brief Advanced static function for the interface 'writeLMem'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void Top_writeLMem_run(
	max_engine_t *engine,
	Top_writeLMem_actions_t *interface_actions);

/**
 * \brief Advanced static non-blocking function for the interface 'writeLMem'.
 *
 * Schedule the actions to run on the engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in] interface_actions Actions to be executed.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Top_writeLMem_run_nonblock(
	max_engine_t *engine,
	Top_writeLMem_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'writeLMem'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void Top_writeLMem_run_group(max_group_t *group, Top_writeLMem_actions_t *interface_actions);

/**
 * \brief Group run advanced static non-blocking function for the interface 'writeLMem'.
 * 
 *
 * Schedule the actions to run on the first device available in the group and return immediately.
 * The status of the run must be checked with ::max_wait. 
 * Note that use of ::max_nowait is prohibited with non-blocking running on groups:
 * see the ::max_run_group_nonblock documentation for more explanation.
 *
 * \param [in] group Group to use.
 * \param [in] interface_actions Actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Top_writeLMem_run_group_nonblock(max_group_t *group, Top_writeLMem_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'writeLMem'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void Top_writeLMem_run_array(max_engarray_t *engarray, Top_writeLMem_actions_t *interface_actions[]);

/**
 * \brief Array run advanced static non-blocking function for the interface 'writeLMem'.
 * 
 *
 * Schedule to run the array of actions on the array of engines, and return immediately.
 * The length of interface_actions must match the size of engarray.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * \param [in] engarray The array of devices to use.
 * \param [in] interface_actions The array of actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Top_writeLMem_run_array_nonblock(max_engarray_t *engarray, Top_writeLMem_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* Top_writeLMem_convert(max_file_t *maxfile, Top_writeLMem_actions_t *interface_actions);



/*----------------------------------------------------------------------------*/
/*---------------------------- Interface default -----------------------------*/
/*----------------------------------------------------------------------------*/



/**
 * \brief Basic static function for the interface 'default'.
 * 
 * \param [in] param_NX Interface Parameter "NX".
 * \param [in] param_inAddr Interface Parameter "inAddr".
 * \param [in] param_inBytes Interface Parameter "inBytes".
 * \param [in] param_inH Interface Parameter "inH".
 * \param [in] param_inW Interface Parameter "inW".
 * \param [in] param_outAddr Interface Parameter "outAddr".
 * \param [in] param_outBytes Interface Parameter "outBytes".
 */
void Top(
	uint32_t param_NX,
	uint32_t param_inAddr,
	uint32_t param_inBytes,
	uint32_t param_inH,
	uint32_t param_inW,
	uint32_t param_outAddr,
	uint32_t param_outBytes);

/**
 * \brief Basic static non-blocking function for the interface 'default'.
 * 
 * Schedule to run on an engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 * 
 * 
 * \param [in] param_NX Interface Parameter "NX".
 * \param [in] param_inAddr Interface Parameter "inAddr".
 * \param [in] param_inBytes Interface Parameter "inBytes".
 * \param [in] param_inH Interface Parameter "inH".
 * \param [in] param_inW Interface Parameter "inW".
 * \param [in] param_outAddr Interface Parameter "outAddr".
 * \param [in] param_outBytes Interface Parameter "outBytes".
 * \return A handle on the execution status, or NULL in case of error.
 */
max_run_t *Top_nonblock(
	uint32_t param_NX,
	uint32_t param_inAddr,
	uint32_t param_inBytes,
	uint32_t param_inH,
	uint32_t param_inW,
	uint32_t param_outAddr,
	uint32_t param_outBytes);

/**
 * \brief Advanced static interface, structure for the engine interface 'default'
 * 
 */
typedef struct { 
	uint32_t param_NX; /**<  [in] Interface Parameter "NX". */
	uint32_t param_inAddr; /**<  [in] Interface Parameter "inAddr". */
	uint32_t param_inBytes; /**<  [in] Interface Parameter "inBytes". */
	uint32_t param_inH; /**<  [in] Interface Parameter "inH". */
	uint32_t param_inW; /**<  [in] Interface Parameter "inW". */
	uint32_t param_outAddr; /**<  [in] Interface Parameter "outAddr". */
	uint32_t param_outBytes; /**<  [in] Interface Parameter "outBytes". */
} Top_actions_t;

/**
 * \brief Advanced static function for the interface 'default'.
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in,out] interface_actions Actions to be executed.
 */
void Top_run(
	max_engine_t *engine,
	Top_actions_t *interface_actions);

/**
 * \brief Advanced static non-blocking function for the interface 'default'.
 *
 * Schedule the actions to run on the engine and return immediately.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * 
 * \param [in] engine The engine on which the actions will be executed.
 * \param [in] interface_actions Actions to be executed.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Top_run_nonblock(
	max_engine_t *engine,
	Top_actions_t *interface_actions);

/**
 * \brief Group run advanced static function for the interface 'default'.
 * 
 * \param [in] group Group to use.
 * \param [in,out] interface_actions Actions to run.
 *
 * Run the actions on the first device available in the group.
 */
void Top_run_group(max_group_t *group, Top_actions_t *interface_actions);

/**
 * \brief Group run advanced static non-blocking function for the interface 'default'.
 * 
 *
 * Schedule the actions to run on the first device available in the group and return immediately.
 * The status of the run must be checked with ::max_wait. 
 * Note that use of ::max_nowait is prohibited with non-blocking running on groups:
 * see the ::max_run_group_nonblock documentation for more explanation.
 *
 * \param [in] group Group to use.
 * \param [in] interface_actions Actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Top_run_group_nonblock(max_group_t *group, Top_actions_t *interface_actions);

/**
 * \brief Array run advanced static function for the interface 'default'.
 * 
 * \param [in] engarray The array of devices to use.
 * \param [in,out] interface_actions The array of actions to run.
 *
 * Run the array of actions on the array of engines.  The length of interface_actions
 * must match the size of engarray.
 */
void Top_run_array(max_engarray_t *engarray, Top_actions_t *interface_actions[]);

/**
 * \brief Array run advanced static non-blocking function for the interface 'default'.
 * 
 *
 * Schedule to run the array of actions on the array of engines, and return immediately.
 * The length of interface_actions must match the size of engarray.
 * The status of the run can be checked either by ::max_wait or ::max_nowait;
 * note that one of these *must* be called, so that associated memory can be released.
 *
 * \param [in] engarray The array of devices to use.
 * \param [in] interface_actions The array of actions to run.
 * \return A handle on the execution status of the actions, or NULL in case of error.
 */
max_run_t *Top_run_array_nonblock(max_engarray_t *engarray, Top_actions_t *interface_actions[]);

/**
 * \brief Converts a static-interface action struct into a dynamic-interface max_actions_t struct.
 *
 * Note that this is an internal utility function used by other functions in the static interface.
 *
 * \param [in] maxfile The maxfile to use.
 * \param [in] interface_actions The interface-specific actions to run.
 * \return The dynamic-interface actions to run, or NULL in case of error.
 */
max_actions_t* Top_convert(max_file_t *maxfile, Top_actions_t *interface_actions);

/**
 * \brief Initialise a maxfile.
 */
max_file_t* Top_init(void);

/* Error handling functions */
int Top_has_errors(void);
const char* Top_get_errors(void);
void Top_clear_errors(void);
/* Free statically allocated maxfile data */
void Top_free(void);
/* returns: -1 = error running command; 0 = no error reported */
int Top_simulator_start(void);
/* returns: -1 = error running command; 0 = no error reported */
int Top_simulator_stop(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* SLIC_DECLARATIONS_Top_H */

